// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/dfm_normal_gamma.h"

#include "core/models/dfm_support.h"
#include "core/models/model_support.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bayests
{

namespace
{

using core::draw_diagonal_precision;
using core::draw_factor_path;
using core::draw_normal_precision;
using core::fill_lagged_factors;
using core::fill_lambda;
using core::identified_loadings;
using core::lambda_row_width;
using core::stacked_response;
using core::transition_residuals;

/// The observed series as the sampler works with them: k x tt, one period per
/// column.
///
/// Built through stacked_response() rather than by transposing `train.y`,
/// because `y` is allowed to arrive already stacked -- a single row or a single
/// column of vec(y') is how the HDF5 files store it -- and only the stacked
/// vector is the same object in all three cases.
arma::mat response_by_period(const TrainData &train, const int k, const int tt)
{
    return arma::reshape(stacked_response(train), k, tt);
}

} // namespace

DfmNormalGammaDraws DfmNormalGammaSampler::draw_coefficients(const DfmNormalGammaInput &input,
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

    const int n_a = input.spec.n_factor_a();
    const bool use_a = n_a > 0;

    DfmNormalGammaDraws out;
    out.lambda = arma::mat(k * n, iterations);
    out.factors = arma::mat(n * tt, iterations);
    out.u_sigma_inv = arma::mat(k, iterations);
    out.v_sigma_inv = arma::mat(n, iterations);

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
        // the serially independent factor this model then has, and its prior on
        // the first state is V. One code path rather than a block-diagonal
        // special case that would draw the same thing.
        a_mat = arma::zeros<arma::mat>(n, n);
    }

    // The order chan_jeliazkov_2009 sees, which is one where this model has
    // none. Everything below that indexes the prior states uses it.
    const int p_state = std::max(p, 1);

    // Error terms. Both are diagonal, so both are carried as the diagonal.
    const arma::vec u_post_shape = input.u_sigma_prior.shape + tt * 0.5;
    const arma::vec &u_prior_rate = input.u_sigma_prior.rate;
    arma::vec u_sigma_inv = input.initial.u_sigma_inv;

    const arma::vec v_post_shape = input.v_sigma_prior.shape + tt * 0.5;
    const arma::vec &v_prior_rate = input.v_sigma_prior.rate;
    arma::vec v_sigma_inv = input.initial.v_sigma_inv;

    arma::mat factors, u, v;

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        // Block 1: Draw the factor path ----
        //
        // Both precisions are carried as their diagonals and the draw wants
        // covariances, so this is where they are turned back -- once per draw,
        // not once per period.
        const arma::mat u_sigma = arma::diagmat(1.0 / u_sigma_inv);
        const arma::mat v_sigma = arma::diagmat(1.0 / v_sigma_inv);

        factors = draw_factor_path(x_t, lambda, u_sigma, v_sigma, a_mat, n, p_state);

        // Block 2: Draw the loadings, equation by equation ----
        //
        // Not as one vector, because the equations do not share a design: row i
        // regresses on the first min(i, N) factors, and while it is inside the
        // identifying block it also carries the fixed unit loading on factor i,
        // which moves to the left-hand side.
        int pos = 0;
        for (int i = 1; i < k; i++)
        {
            const int width = lambda_row_width(i, n);
            const arma::mat f_i = factors.rows(0, width - 1);
            const arma::rowvec response =
                (i < n) ? arma::rowvec(x_t.row(i) - factors.row(i)) : x_t.row(i);

            const arma::mat prior_vinv =
                lambda_prior_vinv.submat(pos, pos, pos + width - 1, pos + width - 1);
            const arma::mat post_v = prior_vinv + (f_i * arma::trans(f_i)) * u_sigma_inv(i);
            const arma::vec rhs = prior_vinv * lambda_prior_mu.subvec(pos, pos + width - 1) +
                                  f_i * arma::trans(response) * u_sigma_inv(i);

            lambda.submat(i, 0, i, width - 1) = arma::trans(draw_normal_precision(post_v, rhs));
            pos += width;
        }

        // Block 3: Draw the idiosyncratic precision ----
        u = x_t - lambda * factors;
        draw_diagonal_precision(u_sigma_inv, u, u_post_shape, u_prior_rate);

        // Block 4: Draw the factor innovation precision ----
        v = transition_residuals(factors, a_mat, n, p);
        draw_diagonal_precision(v_sigma_inv, v, v_post_shape, v_prior_rate);

        // Block 5: Draw the factor transition ----
        //
        // The N equations share their regressors, so the SUR design is
        // kron(X_a', I_N) and its posterior precision collapses to
        // kron(X_a X_a', V^-1) -- the same identity VecKlgs2010 turns on, and
        // the reason no (tt N) x (N^2 p) matrix is built here either.
        if (use_a)
        {
            fill_lagged_factors(x_a, factors, n, p);
            const arma::mat v_prec = arma::diagmat(v_sigma_inv);

            const arma::mat post_v =
                a_prior_vinv + arma::kron(x_a * arma::trans(x_a), v_prec);
            a = draw_normal_precision(
                post_v, a_prior_vinv * a_prior_mu +
                            arma::vectorise(v_prec * (factors * arma::trans(x_a))));
            a_mat = arma::reshape(a, n, n * p);
        }

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            out.lambda.col(draw_pos) = arma::vectorise(lambda);
            out.factors.col(draw_pos) = arma::vectorise(factors);
            out.u_sigma_inv.col(draw_pos) = u_sigma_inv;
            out.v_sigma_inv.col(draw_pos) = v_sigma_inv;
            if (use_a)
            {
                out.a.col(draw_pos) = a;
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws DfmNormalGammaSampler::forecast(const DfmNormalGammaInput &input,
                                              const DfmNormalGammaDraws &coefficients,
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
        const arma::vec u_sd = 1.0 / arma::sqrt(coefficients.u_sigma_inv.col(draw));
        const arma::vec v_sd = 1.0 / arma::sqrt(coefficients.v_sigma_inv.col(draw));

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

arma::mat DfmNormalGammaSampler::log_likelihood(const DfmNormalGammaInput &input,
                                                const DfmNormalGammaDraws &coefficients) const
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

    arma::mat loglik(draws, tt);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        const arma::mat lambda = arma::reshape(coefficients.lambda.col(draw), k, n);
        const arma::mat factors = arma::reshape(coefficients.factors.col(draw), n, tt);
        const arma::vec u_sigma_inv = coefficients.u_sigma_inv.col(draw);

        // U is diagonal, so the determinant term is a sum of logs and the
        // quadratic form is a weighted sum of squares -- no k x k anything.
        const double part_b = arma::accu(arma::log(u_sigma_inv)) / 2;
        const arma::mat u = x_t - lambda * factors;

        for (int i = 0; i < tt; i++)
        {
            const double part_c = -arma::dot(u_sigma_inv, arma::square(u.col(i))) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
