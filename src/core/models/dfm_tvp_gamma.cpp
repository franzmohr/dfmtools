// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/dfm_tvp_gamma.h"

#include "core/algorithms/kalman_durbin_koopman_2002.h"
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
using core::draw_factor_path_tvp;
using core::draw_normal_precision;
using core::fill_lagged_factors;
using core::fill_stacked_loadings;
using core::fill_stacked_transition;
using core::fill_transition_design;
using core::lambda_row_width;
using core::response_by_period;
using core::stacked_identified_loadings;
using core::transition_residuals_tvp;

/// The state variance of a random walk and the state it starts from, drawn
/// together -- the pair every time-varying block here ends with.
///
/// `path` is n x tt and `init` the state before it. The innovations are
/// differences of the path against its own lag, so their sum of squares is what
/// the inverse gamma posterior of the variance adds to the prior rate; `init` is
/// then normal, its only data being the first period of the path.
///
/// `sigma` is the variance itself, not its inverse, and is written in place --
/// which is what the next draw of the path is handed. `post_shape` is the prior
/// shape plus tt/2 and is formed once by the caller.
void draw_random_walk_state(arma::vec &sigma, arma::vec &init, const arma::mat &path,
                            const arma::vec &post_shape, const arma::vec &prior_rate,
                            const NormalPrior &init_prior)
{
    const arma::uword n = path.n_rows;
    const arma::uword tt = path.n_cols;

    arma::mat differences(n, tt);
    differences.col(0) = path.col(0) - init;
    differences.cols(1, tt - 1) = path.cols(1, tt - 1) - path.cols(0, tt - 2);

    const arma::vec sse = arma::sum(arma::square(differences), 1);
    for (arma::uword i = 0; i < n; i++)
    {
        sigma(i) = 1.0 / arma::randg<double>(arma::distr_param(
                             post_shape(i), 1.0 / (prior_rate(i) + sse(i) * 0.5)));
    }

    const arma::mat init_precision = arma::diagmat(1.0 / sigma);
    init = draw_normal_precision(init_prior.v_inv + init_precision,
                                 init_prior.v_inv * init_prior.mu + init_precision * path.col(0));
}

} // namespace

DfmTvpGammaDraws DfmTvpGammaSampler::draw_coefficients(const DfmTvpGammaInput &input,
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
    // the identifying block. Nothing to draw, and nothing to give a state
    // equation to.
    const bool use_lambda = n_lambda > 0;

    DfmTvpGammaDraws out;
    out.lambda = arma::mat(static_cast<arma::uword>(k) * n * tt, iterations);
    out.factors = arma::mat(static_cast<arma::uword>(n) * tt, iterations);
    out.u_sigma_inv = arma::mat(k, iterations);
    out.v_sigma_inv = arma::mat(n, iterations);

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
        // only the cells the Kronecker product is non-zero in, and the rest are
        // structurally zero for the life of the chain.
        z_a = arma::zeros<arma::mat>(static_cast<arma::uword>(n) * tt, n_a);

        out.a = arma::mat(static_cast<arma::uword>(n_a) * tt, iterations);
        out.a_sigma = arma::mat(n_a, iterations);
    }
    else
    {
        // A transition of order zero is a zero transition, which is a thing
        // chan_jeliazkov_2009 can be handed: f_t = 0 f_{t-1} + v_t is exactly
        // the serially independent factor this model then has. One code path
        // rather than a special case, and the same one DfmNormalGamma takes.
        a_stack = arma::zeros<arma::mat>(n, n);
    }

    // The order chan_jeliazkov_2009 sees, which is one where this model has
    // none.
    const int p_state = std::max(p, 1);

    // Both error terms are diagonal, constant, and carried as the diagonal:
    // this model's drift is in the coefficients.
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

        const arma::mat u_sigma = arma::diagmat(1.0 / u_sigma_inv);
        const arma::mat v_sigma = arma::diagmat(1.0 / v_sigma_inv);

        // Block 1: Draw the factor path ----
        //
        // The measurement matrix now differs from period to period, so the
        // assembly's Z' U^-1 Z is formed tt times rather than once. That is what
        // a drifting Lambda costs; see draw_factor_path_tvp().
        factors = draw_factor_path_tvp(x_t, lambda_stack, u_sigma, v_sigma, a_stack, n, p_state);

        if (use_a)
        {
            fill_lagged_factors(x_a, factors, n, p);
        }

        // Block 2: Draw the loading paths, row by row ----
        //
        // Row by row for the reason DfmNormalGamma draws the loadings equation
        // by equation, which the state equation does not change: row i regresses
        // on the first min(i, N) factors and carries min(i, N) free elements, so
        // the rows share no design. Given the factors and a diagonal U they are
        // also conditionally independent, so each row is a state path of its own
        // width rather than a slice of one wide state.
        //
        // While row i is inside the identifying block it carries the fixed unit
        // loading on factor i, which moves to the left-hand side.
        if (use_lambda)
        {
            int pos = 0;
            for (int i = 1; i < k; i++)
            {
                const int width = lambda_row_width(i, n);
                const arma::mat z_i = arma::trans(factors.rows(0, width - 1));
                const arma::mat y_i =
                    (i < n) ? arma::mat(x_t.row(i) - factors.row(i)) : arma::mat(x_t.row(i));

                const arma::mat u_sigma_i(1, 1, arma::fill::value(1.0 / u_sigma_inv(i)));
                const arma::mat sigma_i =
                    arma::diagmat(lambda_sigma.subvec(pos, pos + width - 1));

                lambda_path.rows(pos, pos + width - 1) =
                    kalman_durbin_koopman_2002(y_i, z_i, u_sigma_i, sigma_i,
                                               arma::eye<arma::mat>(width, width),
                                               lambda_init.subvec(pos, pos + width - 1), sigma_i)
                        .cols(0, tt - 1);
                pos += width;
            }

            // Draw the state variance and the loadings before the sample
            draw_random_walk_state(lambda_sigma, lambda_init, lambda_path,
                                   lambda_sigma_post_shape, input.lambda_prior.sigma.rate,
                                   input.lambda_prior.initial_state);

            fill_stacked_loadings(lambda_stack, lambda_path, k, n);
        }

        // Block 3: Draw the idiosyncratic precision ----
        //
        // One measurement matrix per period, so the fitted values are formed
        // period by period rather than as one product.
        u = x_t;
        for (int t = 0; t < tt; t++)
        {
            u.col(t) -= lambda_stack.rows(t * k, (t + 1) * k - 1) * factors.col(t);
        }
        draw_diagonal_precision(u_sigma_inv, u, u_post_shape, u_prior_rate);

        // Block 4: Draw the factor innovation precision ----
        v = use_a ? transition_residuals_tvp(factors, a_stack, x_a, n) : factors;
        draw_diagonal_precision(v_sigma_inv, v, v_post_shape, v_prior_rate);

        // Block 5: Draw the transition path ----
        //
        // As one state of N^2 p elements against the SUR design kron(x_t', I_N).
        // The identity DfmNormalGamma turns on -- the N equations sharing their
        // regressors, so the posterior precision collapses to kron(X X', V^-1)
        // and no (tt N) x (N^2 p) matrix is ever built -- is about a single
        // coefficient vector and does not survive the coefficients becoming a
        // path: each period has its own, and the design has to be spelled out.
        if (use_a)
        {
            fill_transition_design(z_a, x_a, n);

            // Against the precision block 4 has just drawn, not the one the
            // factor path was drawn under at the top of the iteration -- the
            // same conditioning DfmNormalGamma's transition block uses.
            a_path = kalman_durbin_koopman_2002(factors, z_a, arma::diagmat(1.0 / v_sigma_inv),
                                                arma::diagmat(a_sigma), a_B, a_init,
                                                arma::diagmat(a_sigma))
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
            out.u_sigma_inv.col(draw_pos) = u_sigma_inv;
            out.v_sigma_inv.col(draw_pos) = v_sigma_inv;

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

ForecastDraws DfmTvpGammaSampler::forecast(const DfmTvpGammaInput &input,
                                           const DfmTvpGammaDraws &coefficients,
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
    // in-sample period and nothing wider. A caller that hands over the whole
    // path has skipped the slice the io layer makes, and every horizon below
    // would then read the first period's numbers out of it.
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

arma::mat DfmTvpGammaSampler::log_likelihood(const DfmTvpGammaInput &input,
                                             const DfmTvpGammaDraws &coefficients) const
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

    // Every period under its own loadings, so the whole path is wanted here --
    // not the single period a forecast slices out.
    const arma::uword width = static_cast<arma::uword>(k) * n;
    if (coefficients.lambda.n_rows != width * tt)
    {
        throw std::invalid_argument(
            "the pointwise log likelihood of a time-varying model scores every period under its "
            "own loadings, so lambda must have " +
            std::to_string(width * tt) + " rows, got " +
            std::to_string(coefficients.lambda.n_rows));
    }

    arma::mat loglik(draws, tt);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        const arma::mat factors = arma::reshape(coefficients.factors.col(draw), n, tt);
        const arma::vec u_sigma_inv = coefficients.u_sigma_inv.col(draw);

        // U is diagonal and does not move, so the determinant term is a sum of
        // logs formed once per draw and the quadratic form is a weighted sum of
        // squares -- no k x k anything.
        const double part_b = arma::accu(arma::log(u_sigma_inv)) / 2;

        for (int i = 0; i < tt; i++)
        {
            const arma::mat lambda = arma::reshape(
                coefficients.lambda.submat(i * width, draw, (i + 1) * width - 1, draw), k, n);
            const arma::vec u = x_t.col(i) - lambda * factors.col(i);

            const double part_c = -arma::dot(u_sigma_inv, arma::square(u)) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
