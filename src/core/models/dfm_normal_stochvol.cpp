// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/dfm_normal_stochvol.h"

#include "core/algorithms/stochvol_ocsn_2007.h"
#include "core/models/dfm_support.h"
#include "core/models/model_support.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bayests
{

namespace
{

using core::accumulate_transition_moments;
using core::draw_factor_path_sv;
using core::draw_normal_precision;
using core::draw_stochvol_state;
using core::fill_lagged_factors;
using core::fill_lambda;
using core::fill_stacked_diagonal;
using core::identified_loadings;
using core::lambda_row_width;
using core::response_by_period;
using core::transition_residuals;

/// The terminal period of a precision path, whatever the caller brought.
///
/// The forecast holds the volatility at its last in-sample value, so all it ever
/// needs is the last block -- and a host that reads only that block out of a file
/// hands over a matrix `width` rows tall rather than `width * tt`. Counting back
/// from the end covers both without asking which one it was given.
arma::vec terminal_block(const arma::mat &path, const arma::uword draw, const arma::uword width,
                         const char *what)
{
    if (path.n_rows == 0 || path.n_rows % width != 0)
    {
        throw std::invalid_argument(
            std::string("posterior draws of ") + what + " must have a multiple of " +
            std::to_string(width) + " rows, one block per period, got " +
            std::to_string(path.n_rows));
    }
    return path.submat(path.n_rows - width, draw, path.n_rows - 1, draw);
}

} // namespace

DfmNormalStochvolDraws DfmNormalStochvolSampler::draw_coefficients(
    const DfmNormalStochvolInput &input, Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;         // observed series
    const int n = input.spec.n_factors; // factors
    const int p = input.spec.p;         // order of the factor transition
    const int iterations = input.spec.iterations;
    const int burnin = input.spec.burnin;
    const int draws = input.spec.draws();

    const int tt = static_cast<int>(input.train.periods(k));
    const arma::mat x_t = response_by_period(input.train, k, tt);

    const int n_a = input.spec.n_factor_a();
    const bool use_a = n_a > 0;

    DfmNormalStochvolDraws out;
    out.lambda = arma::mat(k * n, iterations);
    out.factors = arma::mat(n * tt, iterations);
    out.u_sigma_inv = arma::mat(k * tt, iterations);
    out.v_sigma_inv = arma::mat(n * tt, iterations);

    // Loadings. The leading N x N block is the identification and is never
    // drawn: ones on the diagonal, zeros above, free below.
    arma::mat lambda = identified_loadings(k, n);
    fill_lambda(lambda, input.initial.lambda);

    const arma::vec &lambda_prior_mu = input.lambda_prior.mu;
    const arma::mat &lambda_prior_vinv = input.lambda_prior.v_inv;

    // Factor transition
    arma::vec a, a_prior_mu;
    arma::mat a_prior_vinv, a_mat, x_a;
    if (use_a)
    {
        a_prior_mu = input.a_prior.mu;
        a_prior_vinv = input.a_prior.v_inv;
        a = input.initial.a;
        a_mat = arma::reshape(a, n, n * p);
        x_a = arma::mat(n * p, tt);
        out.a = arma::mat(n_a, iterations);
    }
    else
    {
        // A transition of order zero is a zero transition, which is a thing
        // chan_jeliazkov_2009 can be handed: f_t = 0 f_{t-1} + v_t is exactly
        // the serially independent -- but still heteroskedastic -- factor this
        // model then has. One code path rather than a special case.
        a_mat = arma::zeros<arma::mat>(n, n);
    }

    // The order chan_jeliazkov_2009 sees, which is one where this model has
    // none. Everything below that indexes the prior states uses it.
    const int p_state = std::max(p, 1);

    // Idiosyncratic volatility
    arma::mat u_h = input.initial.u_h; // tt x k
    arma::vec u_h_init = input.initial.u_h_init;
    arma::vec u_h_sigma = input.initial.u_h_sigma;
    const arma::vec &u_h_offset = input.u_sigma_prior.offset;
    const arma::vec u_h_sigma_post_shape = input.u_sigma_prior.state.sigma.shape + tt * 0.5;
    const arma::vec &u_h_sigma_prior_rate = input.u_sigma_prior.state.sigma.rate;

    // Factor innovation volatility
    arma::mat v_h = input.initial.v_h; // tt x n
    arma::vec v_h_init = input.initial.v_h_init;
    arma::vec v_h_sigma = input.initial.v_h_sigma;
    const arma::vec &v_h_offset = input.v_sigma_prior.offset;
    const arma::vec v_h_sigma_post_shape = input.v_sigma_prior.state.sigma.shape + tt * 0.5;
    const arma::vec &v_h_sigma_prior_rate = input.v_sigma_prior.state.sigma.rate;

    // The two covariance stacks chan_jeliazkov_2009 reads, one K x K or N x N
    // block per period. Zeroed once here and never zeroed again: both
    // covariances are diagonal, so fill_stacked_diagonal only ever writes the
    // diagonal and everything else stays at the zero it starts with. Allocated
    // outside the loop as well -- at 100 series over 300 periods this is 24 MB,
    // which is worth allocating once and worth not touching three million cells
    // of per draw.
    arma::mat u_stack(k * tt, k, arma::fill::zeros);
    arma::mat v_stack(n * tt, n, arma::fill::zeros);

    arma::mat u_precision = arma::exp(-u_h); // tt x k
    arma::mat v_precision = arma::exp(-v_h); // tt x n
    fill_stacked_diagonal(u_stack, arma::exp(u_h));
    fill_stacked_diagonal(v_stack, arma::exp(v_h));

    arma::mat factors, u, v;

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        // Block 1: Draw the factor path ----
        //
        // Unlike DfmNormalGamma there is no per-draw inversion here: the stacks
        // hold covariances already, and chan_jeliazkov_2009 takes their
        // precisions itself -- one diagonal reciprocal per period, since it
        // recognises a diagonal covariance and does not factorise one.
        factors = draw_factor_path_sv(x_t, lambda, u_stack, v_stack, a_mat, n, p_state);

        // Block 2: Draw the loadings, equation by equation ----
        //
        // Not as one vector, because the equations do not share a design: row i
        // regresses on the first min(i, N) factors, and while it is inside the
        // identifying block it also carries the fixed unit loading on factor i,
        // which moves to the left-hand side.
        //
        // Weighted rather than scaled, which is the whole point of the model.
        // Where DfmNormalGamma multiplies the cross-products by the single
        // number u_sigma_inv(i), row i here carries one weight per period, so
        // the periods in which series i was quiet identify its loadings and the
        // periods in which it was wild largely do not.
        int pos = 0;
        for (int i = 1; i < k; i++)
        {
            const int width = lambda_row_width(i, n);
            const arma::mat f_i = factors.rows(0, width - 1);
            const arma::rowvec response =
                (i < n) ? arma::rowvec(x_t.row(i) - factors.row(i)) : arma::rowvec(x_t.row(i));

            arma::mat f_weighted = f_i;
            f_weighted.each_row() %= arma::trans(u_precision.col(i));

            const arma::mat prior_vinv =
                lambda_prior_vinv.submat(pos, pos, pos + width - 1, pos + width - 1);
            const arma::mat post_v = prior_vinv + f_weighted * arma::trans(f_i);
            const arma::vec rhs = prior_vinv * lambda_prior_mu.subvec(pos, pos + width - 1) +
                                  f_weighted * arma::trans(response);

            lambda.submat(i, 0, i, width - 1) = arma::trans(draw_normal_precision(post_v, rhs));
            pos += width;
        }

        // Block 3: Draw the idiosyncratic log-volatility ----
        //
        // The factored routine, which is the one test/unit_stochvol.cpp covers
        // and which draws each path with a banded Cholesky rather than
        // factorising a dense tt x tt precision. The k series are independent
        // given the factors, so one call handles all of them.
        u = x_t - lambda * factors;
        u_h = stochvol_ocsn_2007(arma::trans(u), u_h, u_h_sigma, u_h_init, u_h_offset);
        draw_stochvol_state(u_h_sigma, u_h_init, u_h, u_h_sigma_post_shape, u_h_sigma_prior_rate,
                            input.u_sigma_prior.state.initial_state);

        const arma::mat u_variance = arma::exp(u_h);
        u_precision = 1.0 / u_variance;
        fill_stacked_diagonal(u_stack, u_variance);

        // Block 4: Draw the factor innovation log-volatility ----
        //
        // The residual of period t is the transition's, including the first p
        // periods, where it is the transition truncated at the zero factors
        // before the sample. Those are as much a draw from N(0, V_t) as the rest,
        // so all tt of them inform the volatility.
        v = transition_residuals(factors, a_mat, n, p);
        v_h = stochvol_ocsn_2007(arma::trans(v), v_h, v_h_sigma, v_h_init, v_h_offset);
        draw_stochvol_state(v_h_sigma, v_h_init, v_h, v_h_sigma_post_shape, v_h_sigma_prior_rate,
                            input.v_sigma_prior.state.initial_state);

        const arma::mat v_variance = arma::exp(v_h);
        v_precision = 1.0 / v_variance;
        fill_stacked_diagonal(v_stack, v_variance);

        // Block 5: Draw the factor transition ----
        //
        // DfmNormalGamma's Kronecker collapse is gone. Its posterior precision
        // is kron(X X', V^-1) because V is the same in every period, and
        // sum_t kron(x_t x_t', V_t^-1) does not factor that way.
        //
        // What survives is that V_t is diagonal, which is what
        // accumulate_transition_moments turns on -- see there for the identity
        // and for the scatter. The prior, which need not respect that structure,
        // is added whole. No (tt N) x (N^2 p) matrix is built here either.
        if (use_a)
        {
            fill_lagged_factors(x_a, factors, n, p);

            arma::mat a_post_v = a_prior_vinv;
            arma::vec a_rhs = a_prior_vinv * a_prior_mu;
            accumulate_transition_moments(a_post_v, a_rhs, x_a, factors, v_precision, n, p);

            a = draw_normal_precision(a_post_v, a_rhs);
            a_mat = arma::reshape(a, n, n * p);
        }

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            out.lambda.col(draw_pos) = arma::vectorise(lambda);
            out.factors.col(draw_pos) = arma::vectorise(factors);

            // Periods stacked within a column, k or n values per period, which
            // is the transpose of the tt x k the volatility block works in.
            out.u_sigma_inv.col(draw_pos) = arma::vectorise(arma::trans(u_precision));
            out.v_sigma_inv.col(draw_pos) = arma::vectorise(arma::trans(v_precision));

            if (use_a)
            {
                out.a.col(draw_pos) = a;
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws DfmNormalStochvolSampler::forecast(const DfmNormalStochvolInput &input,
                                                 const DfmNormalStochvolDraws &coefficients,
                                                 Reporter &reporter) const
{
    const int k = input.spec.k;
    const int n = input.spec.n_factors;
    const int p = input.spec.p;
    const int h = input.spec.h;
    const bool use_a = input.use_a();

    if (k <= 0 || n <= 0)
    {
        throw std::invalid_argument("model must have at least one observed series (k) and one "
                                    "factor (n_factors)");
    }
    if (h <= 0)
    {
        throw std::invalid_argument("forecast horizon (h) must be positive");
    }
    if (!coefficients.has_factors())
    {
        throw std::invalid_argument("posterior draws of the factors are missing; a dynamic factor "
                                    "model forecasts by running the transition on from them");
    }
    if (coefficients.lambda.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of lambda are missing");
    }
    if (coefficients.u_sigma_inv.n_elem == 0 || coefficients.v_sigma_inv.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of the error precisions are missing");
    }
    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("the factors have a transition of order " +
                                    std::to_string(p) + " but posterior draws of a are missing");
    }

    const int tt = static_cast<int>(input.train.periods(k));
    if (tt < p)
    {
        throw std::invalid_argument("a transition of order " + std::to_string(p) +
                                    " needs that many factors to start from, and the sample has " +
                                    std::to_string(tt));
    }
    if (coefficients.factors.n_rows != static_cast<arma::uword>(n * tt))
    {
        throw std::invalid_argument(
            "posterior draws of the factors must have " + std::to_string(n * tt) +
            " rows for " + std::to_string(n) + " factors over " + std::to_string(tt) +
            " periods, got " + std::to_string(coefficients.factors.n_rows));
    }

    const arma::uword draws = coefficients.iterations();
    arma::mat fcst(h * k, draws);

    // The path carries the p factors the transition needs before the first
    // horizon, so column p + i is horizon i and the lag lookup is one expression
    // at every horizon rather than a split between what is history and what is
    // already forecast.
    arma::mat path(n, p + h);

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(draws));

        const arma::mat lambda = arma::reshape(coefficients.lambda.col(draw), k, n);

        // Both volatilities held at the last in-sample period. See the header
        // for why that is the convention rather than a simulation forward.
        const arma::vec u_sd = 1.0 / arma::sqrt(terminal_block(
                                        coefficients.u_sigma_inv, draw,
                                        static_cast<arma::uword>(k), "u_sigma_inv"));
        const arma::vec v_sd = 1.0 / arma::sqrt(terminal_block(
                                        coefficients.v_sigma_inv, draw,
                                        static_cast<arma::uword>(n), "v_sigma_inv"));

        if (p > 0)
        {
            const arma::mat drawn = arma::reshape(coefficients.factors.col(draw), n, tt);
            path.head_cols(p) = drawn.tail_cols(p);
        }

        const arma::mat a_mat =
            use_a ? arma::reshape(coefficients.a.col(draw), n, n * p) : arma::mat();

        for (int i = 0; i < h; i++)
        {
            arma::vec f = v_sd % arma::randn<arma::vec>(n);
            for (int j = 1; j <= p; j++)
            {
                f += a_mat.cols((j - 1) * n, j * n - 1) * path.col(p + i - j);
            }
            path.col(p + i) = f;

            fcst.submat(i * k, draw, (i + 1) * k - 1, draw) =
                lambda * f + u_sd % arma::randn<arma::vec>(k);
        }
    }

    reporter.finish();
    return ForecastDraws{fcst};
}

arma::mat DfmNormalStochvolSampler::log_likelihood(const DfmNormalStochvolInput &input,
                                                   const DfmNormalStochvolDraws &coefficients) const
{
    const int k = input.spec.k;
    const int n = input.spec.n_factors;

    if (k <= 0 || n <= 0)
    {
        throw std::invalid_argument("model must have at least one observed series (k) and one "
                                    "factor (n_factors)");
    }
    if (coefficients.u_sigma_inv.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of u_sigma_inv are missing");
    }
    if (!coefficients.has_factors())
    {
        throw std::invalid_argument("posterior draws of the factors are missing; the likelihood of "
                                    "a dynamic factor model is conditional on them");
    }
    if (coefficients.lambda.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of lambda are missing");
    }

    const int tt = static_cast<int>(input.train.periods(k));
    const arma::mat x_t = response_by_period(input.train, k, tt);
    const arma::uword draws = coefficients.iterations();

    // Unlike the forecast, this needs the whole path: every period has its own
    // precision and every period contributes a column.
    if (coefficients.u_sigma_inv.n_rows != static_cast<arma::uword>(k * tt))
    {
        throw std::invalid_argument(
            "posterior draws of u_sigma_inv must have " + std::to_string(k * tt) +
            " rows for " + std::to_string(k) + " series over " + std::to_string(tt) +
            " periods, got " + std::to_string(coefficients.u_sigma_inv.n_rows));
    }

    arma::mat loglik(draws, tt);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        const arma::mat lambda = arma::reshape(coefficients.lambda.col(draw), k, n);
        const arma::mat factors = arma::reshape(coefficients.factors.col(draw), n, tt);
        const arma::mat u = x_t - lambda * factors;

        // U_t is diagonal, so the determinant term is a sum of logs and the
        // quadratic form is a weighted sum of squares -- no k x k anything. The
        // determinant does move into the loop, though, which is the one thing
        // this does not share with DfmNormalGamma: there it is one number per
        // draw, here it is one per period.
        const arma::mat u_sigma_inv =
            arma::reshape(coefficients.u_sigma_inv.col(draw), k, tt);

        for (int i = 0; i < tt; i++)
        {
            const arma::vec precision = u_sigma_inv.col(i);
            const double part_b = arma::accu(arma::log(precision)) / 2;
            const double part_c = -arma::dot(precision, arma::square(u.col(i))) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
