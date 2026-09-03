// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/dfm_tvp_stochvol.h"

#include "core/algorithms/kalman_durbin_koopman_2002.h"
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

using core::draw_factor_path;
using core::draw_random_walk_state;
using core::draw_stochvol_state;
using core::fill_lagged_factors;
using core::fill_stacked_diagonal;
using core::fill_stacked_loadings;
using core::fill_stacked_transition;
using core::fill_transition_design;
using core::lambda_row_width;
using core::response_by_period;
using core::stacked_identified_loadings;
using core::transition_residuals_tvp;

/// The terminal period of a precision path, whatever the caller brought.
///
/// The forecast holds the volatility at its last in-sample value, so all it ever
/// needs is the last block -- and a host that reads only that block out of a file
/// hands over a matrix `width` rows tall rather than `width * tt`. Counting back
/// from the end covers both without asking which one it was given. The same
/// helper DfmNormalStochvol carries, for the same pair of shapes.
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

DfmTvpStochvolDraws DfmTvpStochvolSampler::draw_coefficients(const DfmTvpStochvolInput &input,
                                                             Reporter &reporter) const
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

    const int n_lambda = input.spec.n_lambda();
    const int n_a = input.spec.n_factor_a();
    const bool use_a = n_a > 0;

    // A single observed series has no free loading: the whole of Lambda is then
    // the identifying block.
    const bool use_lambda = n_lambda > 0;

    DfmTvpStochvolDraws out;
    out.lambda = arma::mat(static_cast<arma::uword>(k) * n * tt, iterations);
    out.factors = arma::mat(static_cast<arma::uword>(n) * tt, iterations);
    out.u_sigma_inv = arma::mat(static_cast<arma::uword>(k) * tt, iterations);
    out.v_sigma_inv = arma::mat(static_cast<arma::uword>(n) * tt, iterations);

    // The loadings: a path per free element, and the identifying block that is
    // the same in every period and is never drawn.
    arma::mat lambda_stack = stacked_identified_loadings(k, n, tt);
    arma::mat lambda_path;
    arma::vec lambda_sigma, lambda_init, lambda_sigma_post_shape;
    if (use_lambda)
    {
        lambda_path = input.initial.lambda;
        lambda_sigma = 1.0 / input.initial.lambda_sigma_inv.diag();
        lambda_init = input.initial.lambda_init;
        lambda_sigma_post_shape = input.lambda_prior.sigma.shape + tt * 0.5;
        out.lambda_sigma = arma::mat(n_lambda, iterations);
        fill_stacked_loadings(lambda_stack, lambda_path, k, n);
    }

    // The factor transition: a path of vec([A_1 .. A_p]), the SUR design it is
    // drawn against, and the lagged factors that design is built from.
    arma::mat a_path, a_stack, a_B, x_a, z_a;
    arma::vec a_sigma, a_init, a_sigma_post_shape;
    if (use_a)
    {
        a_path = input.initial.a;
        a_sigma = 1.0 / input.initial.a_sigma_inv.diag();
        a_init = input.initial.a_init;
        a_sigma_post_shape = input.a_prior.sigma.shape + tt * 0.5;
        a_B = arma::eye<arma::mat>(n_a, n_a);

        a_stack = arma::mat(static_cast<arma::uword>(n) * tt, n * p);
        fill_stacked_transition(a_stack, a_path, n, p);

        x_a = arma::mat(n * p, tt);

        // Zeroed once and never zeroed again: fill_transition_design() writes
        // only the cells the Kronecker product is non-zero in.
        z_a = arma::zeros<arma::mat>(static_cast<arma::uword>(n) * tt, n_a);

        out.a = arma::mat(static_cast<arma::uword>(n_a) * tt, iterations);
        out.a_sigma = arma::mat(n_a, iterations);
    }
    else
    {
        // A transition of order zero is a zero transition, which is a thing
        // chan_jeliazkov_2009 can be handed: f_t = 0 f_{t-1} + v_t is exactly
        // the serially independent -- but still heteroskedastic -- factor this
        // model then has.
        a_stack = arma::zeros<arma::mat>(n, n);
    }

    // The order chan_jeliazkov_2009 sees, which is one where this model has
    // none.
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

    // The two covariance stacks the factor path reads, one K x K or N x N block
    // per period. Zeroed once and never zeroed again: both covariances are
    // diagonal, so fill_stacked_diagonal only ever writes the diagonal.
    arma::mat u_stack(static_cast<arma::uword>(k) * tt, k, arma::fill::zeros);
    arma::mat v_stack(static_cast<arma::uword>(n) * tt, n, arma::fill::zeros);

    arma::mat u_variance = arma::exp(u_h); // tt x k
    arma::mat u_precision = 1.0 / u_variance;
    arma::mat v_variance = arma::exp(v_h); // tt x n
    arma::mat v_precision = 1.0 / v_variance;
    fill_stacked_diagonal(u_stack, u_variance);
    fill_stacked_diagonal(v_stack, v_variance);

    arma::mat factors, u, v;

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        // Block 1: Draw the factor path ----
        //
        // The one call in this library that stacks all four of the band
        // sampler's per-period arguments: a loading matrix, a measurement
        // covariance, a transition and a transition covariance, each its own per
        // period. draw_factor_path() carries the shift the last two need and the
        // reason the first two do not.
        factors = draw_factor_path(x_t, lambda_stack, u_stack, v_stack, a_stack, n, p_state);

        if (use_a)
        {
            fill_lagged_factors(x_a, factors, n, p);
        }

        // Block 2: Draw the loading paths, row by row ----
        //
        // Row by row for the reason DfmNormalGamma draws the loadings equation
        // by equation: row i regresses on the first min(i, N) factors, so the
        // rows share no design, and given the factors and a diagonal U they are
        // conditionally independent.
        //
        // Weighted per period as well as drifting, which is what neither of the
        // two models this one sits between has. Row i carries one weight per
        // period, so the periods in which series i was quiet identify its
        // loading path and the periods in which it was wild largely do not --
        // and the path is free to move between them, which is exactly the
        // confusion the two halves of this model exist to tell apart.
        if (use_lambda)
        {
            int pos = 0;
            for (int i = 1; i < k; i++)
            {
                const int width = lambda_row_width(i, n);
                const arma::mat z_i = arma::trans(factors.rows(0, width - 1));
                const arma::mat y_i =
                    (i < n) ? arma::mat(x_t.row(i) - factors.row(i)) : arma::mat(x_t.row(i));

                // One 1 x 1 covariance per period, which is the stacked shape
                // the smoother reads for a single-row measurement.
                const arma::mat u_sigma_i = u_variance.col(i);
                const arma::mat sigma_i =
                    arma::diagmat(lambda_sigma.subvec(pos, pos + width - 1));

                lambda_path.rows(pos, pos + width - 1) =
                    kalman_durbin_koopman_2002(y_i, z_i, u_sigma_i, sigma_i,
                                               arma::eye<arma::mat>(width, width),
                                               lambda_init.subvec(pos, pos + width - 1), sigma_i)
                        .cols(0, tt - 1);
                pos += width;
            }

            draw_random_walk_state(lambda_sigma, lambda_init, lambda_path,
                                   lambda_sigma_post_shape, input.lambda_prior.sigma.rate,
                                   input.lambda_prior.initial_state);

            fill_stacked_loadings(lambda_stack, lambda_path, k, n);
        }

        // Block 3: Draw the idiosyncratic log-volatility ----
        //
        // The fitted values are formed period by period, the measurement matrix
        // being one per period. Everything after that is DfmNormalStochvol's
        // block unchanged: the k series are independent given the factors, so one
        // call to the factored routine handles all of them.
        u = x_t;
        for (int t = 0; t < tt; t++)
        {
            u.col(t) -= lambda_stack.rows(t * k, (t + 1) * k - 1) * factors.col(t);
        }
        u_h = stochvol_ocsn_2007(arma::trans(u), u_h, u_h_sigma, u_h_init, u_h_offset);
        draw_stochvol_state(u_h_sigma, u_h_init, u_h, u_h_sigma_post_shape, u_h_sigma_prior_rate,
                            input.u_sigma_prior.state.initial_state);

        u_variance = arma::exp(u_h);
        u_precision = 1.0 / u_variance;
        fill_stacked_diagonal(u_stack, u_variance);

        // Block 4: Draw the factor innovation log-volatility ----
        //
        // The residual of period t is the transition's, including the first p
        // periods, where it is the transition truncated at the zero factors
        // before the sample. Those are as much a draw from N(0, V_t) as the rest,
        // so all tt of them inform the volatility.
        v = use_a ? transition_residuals_tvp(factors, a_stack, x_a, n) : factors;
        v_h = stochvol_ocsn_2007(arma::trans(v), v_h, v_h_sigma, v_h_init, v_h_offset);
        draw_stochvol_state(v_h_sigma, v_h_init, v_h, v_h_sigma_post_shape, v_h_sigma_prior_rate,
                            input.v_sigma_prior.state.initial_state);

        v_variance = arma::exp(v_h);
        v_precision = 1.0 / v_variance;
        fill_stacked_diagonal(v_stack, v_variance);

        // Block 5: Draw the transition path ----
        //
        // As one state of N^2 p elements against the SUR design kron(x_t', I_N),
        // measured with the covariance block 4 has just drawn -- so the smoother
        // sees a heteroskedastic measurement equation, which is the stacked shape
        // v_stack already is.
        //
        // Neither of the two economies the constant-coefficient models take
        // survives here. DfmNormalGamma's Kronecker collapse is a statement about
        // a single coefficient vector, and DfmNormalStochvol's
        // accumulate_transition_moments() is about a single vector under a moving
        // covariance; a path is neither, so the design is spelled out.
        if (use_a)
        {
            fill_transition_design(z_a, x_a, n);

            a_path = kalman_durbin_koopman_2002(factors, z_a, v_stack, arma::diagmat(a_sigma),
                                                a_B, a_init, arma::diagmat(a_sigma))
                         .cols(0, tt - 1);

            draw_random_walk_state(a_sigma, a_init, a_path, a_sigma_post_shape,
                                   input.a_prior.sigma.rate, input.a_prior.initial_state);

            fill_stacked_transition(a_stack, a_path, n, p);
        }

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            // Lambda, one vectorised M x N block per period. Not one vectorise()
            // of the stack: that is (M tt) x N and column-major, so it would
            // interleave the periods rather than stack them.
            const arma::uword width = static_cast<arma::uword>(k) * n;
            for (int t = 0; t < tt; t++)
            {
                out.lambda.submat(t * width, draw_pos, (t + 1) * width - 1, draw_pos) =
                    arma::vectorise(lambda_stack.rows(t * k, (t + 1) * k - 1));
            }

            out.factors.col(draw_pos) = arma::vectorise(factors);

            // Periods stacked within a column, k or n values per period, which
            // is the transpose of the tt x k the volatility block works in.
            out.u_sigma_inv.col(draw_pos) = arma::vectorise(arma::trans(u_precision));
            out.v_sigma_inv.col(draw_pos) = arma::vectorise(arma::trans(v_precision));

            if (use_lambda)
            {
                out.lambda_sigma.col(draw_pos) = lambda_sigma;
            }
            if (use_a)
            {
                out.a.col(draw_pos) = arma::vectorise(a_path);
                out.a_sigma.col(draw_pos) = a_sigma;
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws DfmTvpStochvolSampler::forecast(const DfmTvpStochvolInput &input,
                                              const DfmTvpStochvolDraws &coefficients,
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

    // The coefficients move, so what a forecast starts from is their last
    // in-sample period and nothing wider. The two precisions are read the other
    // way -- terminal_block() takes the last block of whatever it is handed --
    // because a path of them is what the log likelihood also wants, and cutting
    // one to a period is unambiguous where cutting the other would not be.
    if (coefficients.lambda.n_rows != static_cast<arma::uword>(k) * n)
    {
        throw std::invalid_argument(
            "forecasting a time-varying model starts from the last in-sample period, so lambda "
            "must have " +
            std::to_string(k * n) + " rows, got " + std::to_string(coefficients.lambda.n_rows));
    }
    if (use_a && coefficients.a.n_rows != static_cast<arma::uword>(input.spec.n_factor_a()))
    {
        throw std::invalid_argument(
            "forecasting a time-varying model starts from the last in-sample period, so a must "
            "have " +
            std::to_string(input.spec.n_factor_a()) + " rows, got " +
            std::to_string(coefficients.a.n_rows));
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
            "posterior draws of the factors must have " + std::to_string(n * tt) + " rows for " +
            std::to_string(n) + " factors over " + std::to_string(tt) + " periods, got " +
            std::to_string(coefficients.factors.n_rows));
    }

    const arma::uword draws = coefficients.iterations();
    arma::mat fcst(h * k, draws);

    // The path carries the p factors the transition needs before the first
    // horizon, so column p + i is horizon i.
    arma::mat path(n, p + h);

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(draws));

        const arma::mat lambda = arma::reshape(coefficients.lambda.col(draw), k, n);

        // Both volatilities held at the last in-sample period, as every
        // stochastic volatility model here does: the variance of the
        // log-volatility innovations is a state of the chain rather than
        // something the draws carry, so there is nothing to extrapolate the
        // random walk with.
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

arma::mat DfmTvpStochvolSampler::log_likelihood(const DfmTvpStochvolInput &input,
                                                const DfmTvpStochvolDraws &coefficients) const
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

    // Both paths whole: every period is scored under its own loadings and its own
    // precision, and neither the terminal slice a forecast takes nor a single
    // matrix would say what this model claims.
    const arma::uword width = static_cast<arma::uword>(k) * n;
    if (coefficients.lambda.n_rows != width * tt)
    {
        throw std::invalid_argument(
            "the pointwise log likelihood of a time-varying model scores every period under its "
            "own loadings, so lambda must have " +
            std::to_string(width * tt) + " rows, got " +
            std::to_string(coefficients.lambda.n_rows));
    }
    if (coefficients.u_sigma_inv.n_rows != static_cast<arma::uword>(k * tt))
    {
        throw std::invalid_argument(
            "posterior draws of u_sigma_inv must have " + std::to_string(k * tt) + " rows for " +
            std::to_string(k) + " series over " + std::to_string(tt) + " periods, got " +
            std::to_string(coefficients.u_sigma_inv.n_rows));
    }

    arma::mat loglik(draws, tt);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        const arma::mat factors = arma::reshape(coefficients.factors.col(draw), n, tt);
        const arma::mat u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.col(draw), k, tt);

        // U_t is diagonal, so the determinant term is a sum of logs and the
        // quadratic form a weighted sum of squares -- no k x k anything. Both
        // move inside the period loop here, which is what makes this the one of
        // the four that shares neither model's shortcut: the fitted value changes
        // with Lambda_t and the determinant with U_t.
        for (int i = 0; i < tt; i++)
        {
            const arma::mat lambda = arma::reshape(
                coefficients.lambda.submat(i * width, draw, (i + 1) * width - 1, draw), k, n);
            const arma::vec u = x_t.col(i) - lambda * factors.col(i);
            const arma::vec precision = u_sigma_inv.col(i);

            const double part_b = arma::accu(arma::log(precision)) / 2;
            const double part_c = -arma::dot(precision, arma::square(u)) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
