// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/favar_normal_wishart.h"

#include "core/algorithms/wishart.h"
#include "core/models/favar_support.h"
#include "core/models/model_support.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bayests
{

namespace
{

using core::draw_conditional_factor_path;
using core::draw_diagonal_precision;
using core::draw_normal_precision;
using core::favar_identified_loadings;
using core::favar_lambda_row_width;
using core::fill_favar_lambda;
using core::fill_lagged_factors;
using core::obs_factors_by_period;
using core::response_by_period;
using core::stacked_state;
using core::transition_residuals;

/// The symmetric square root of a covariance, for drawing an innovation from it.
///
/// The eigen route rather than a Cholesky, which is the spelling
/// `var_normal_wishart.cpp` uses for the same job. Taken once per draw here
/// rather than once per horizon: Q does not move over the forecast, so the
/// factorisation does not either.
arma::mat covariance_root(const arma::mat &precision)
{
    arma::vec eigval;
    arma::mat eigvec;
    const arma::mat covariance =
        arma::solve(precision, arma::eye<arma::mat>(precision.n_rows, precision.n_rows));
    if (!arma::eig_sym(eigval, eigvec, arma::symmatu(covariance)))
    {
        throw std::runtime_error("the drawn state innovation precision has no symmetric square "
                                 "root; the chain has degenerated");
    }
    return eigvec * arma::diagmat(arma::sqrt(arma::abs(eigval))) * arma::trans(eigvec);
}

} // namespace

FavarNormalWishartDraws
FavarNormalWishartSampler::draw_coefficients(const FavarNormalWishartInput &input,
                                             Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;               // observed series in the panel
    const int n = input.spec.n_factors;       // unobserved factors
    const int n_obs = input.spec.n_obs_factors; // observed factors
    const int ns = input.spec.n_state();      // the state the transition runs over
    const int p = input.spec.p;               // order of the state transition
    const int iterations = input.spec.iterations;
    const int burnin = input.spec.burnin;
    const int draws = input.spec.draws();

    const int tt = static_cast<int>(input.train.periods(k));
    const arma::mat x_t = response_by_period(input.train, k, tt);

    // The observed half of the state. Transposed once here rather than once per
    // draw: it is data and does not change over the chain.
    const arma::mat obs = obs_factors_by_period(input.train);

    const int n_a = input.spec.n_favar_a();
    const bool use_a = n_a > 0;

    FavarNormalWishartDraws out;
    out.lambda = arma::mat(k * ns, iterations);
    out.factors = arma::mat(n * tt, iterations);
    out.u_sigma_inv = arma::mat(k, iterations);
    out.v_sigma_inv = arma::mat(ns * ns, iterations);

    // Loadings. The leading n x n block of the factor columns is unit lower
    // triangular and the observed columns of those rows are zero; neither is
    // ever drawn.
    arma::mat lambda = favar_identified_loadings(k, n, n_obs);
    fill_favar_lambda(lambda, input.initial.lambda, n, n_obs);

    const arma::vec &lambda_prior_mu = input.lambda_prior.mu;
    const arma::mat &lambda_prior_vinv = input.lambda_prior.v_inv;

    // State transition, a VAR over the whole state rather than over the factors
    // alone -- which is the model.
    arma::vec a, a_prior_mu;
    arma::mat a_prior_vinv, a_mat, x_a;
    if (use_a)
    {
        a_prior_mu = input.a_prior.mu;
        a_prior_vinv = input.a_prior.v_inv;
        a = input.initial.a;
        a_mat = arma::reshape(a, ns, ns * p);
        x_a = arma::mat(ns * p, tt);
        out.a = arma::mat(n_a, iterations);
    }
    else
    {
        // A transition of order zero is a zero transition, exactly as in a
        // dynamic factor model: s_t = 0 s_{t-1} + v_t is the serially
        // independent state the model then has, and its prior on the first
        // state is Q.
        a_mat = arma::zeros<arma::mat>(ns, ns);
    }

    // The order chan_jeliazkov_2009 sees, which is one where this model has
    // none.
    const int p_state = std::max(p, 1);

    // Idiosyncratic precision. Diagonal by assumption -- that assumption is what
    // makes this a factor model -- so it is carried as the diagonal.
    const arma::vec u_post_shape = input.u_sigma_prior.shape + tt * 0.5;
    const arma::vec &u_prior_rate = input.u_sigma_prior.rate;
    arma::vec u_sigma_inv = input.initial.u_sigma_inv;

    // State innovation precision. A whole matrix, which is what separates this
    // model from every DFM: see the class comment.
    const arma::mat &v_prior_scale = input.v_sigma_prior.scale;
    const int v_post_df = input.v_sigma_prior.df + tt;
    arma::mat v_sigma_inv = input.initial.v_sigma_inv;

    const arma::mat state_identity = arma::eye<arma::mat>(ns, ns);
    arma::mat factors, state, u, v;

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        // Block 1: Draw the path of the unobserved factors ----
        //
        // The observed factors are handed over as `obs` and held there. They are
        // not subtracted out of the panel first: their contribution to the
        // measurement is one of the cross terms the conditioning removes, so
        // doing it by hand as well would subtract it twice.
        const arma::mat u_sigma = arma::diagmat(1.0 / u_sigma_inv);
        const arma::mat v_sigma = arma::solve(v_sigma_inv, state_identity);

        factors = draw_conditional_factor_path(x_t, lambda, u_sigma, arma::symmatu(v_sigma),
                                               a_mat, obs, ns, p_state);
        state = stacked_state(factors, obs);

        // Block 2: Draw the loadings, equation by equation ----
        //
        // Not as one vector, because the equations do not share a design: the
        // first n rows are the identifying block and carry no free loading at
        // all, and every row after them regresses on the whole state. Unlike a
        // dynamic factor model there is no partial row and nothing moves to the
        // left-hand side -- the identifying block being the identity rather than
        // a unit triangle is exactly what removes that case.
        int pos = 0;
        for (int i = 1; i < k; i++)
        {
            const int width = favar_lambda_row_width(i, n, n_obs);
            if (width == 0)
            {
                continue;
            }

            const arma::mat s_i = state.rows(0, width - 1);
            const arma::rowvec response = x_t.row(i);

            const arma::mat prior_vinv =
                lambda_prior_vinv.submat(pos, pos, pos + width - 1, pos + width - 1);
            const arma::mat post_v = prior_vinv + (s_i * arma::trans(s_i)) * u_sigma_inv(i);
            const arma::vec rhs = prior_vinv * lambda_prior_mu.subvec(pos, pos + width - 1) +
                                  s_i * arma::trans(response) * u_sigma_inv(i);

            lambda.submat(i, 0, i, width - 1) = arma::trans(draw_normal_precision(post_v, rhs));
            pos += width;
        }

        // Block 3: Draw the idiosyncratic precision ----
        u = x_t - lambda * state;
        draw_diagonal_precision(u_sigma_inv, u, u_post_shape, u_prior_rate);

        // Block 4: Draw the state innovation precision ----
        //
        // The one block that is a matrix draw rather than a diagonal one. Q
        // couples the factor innovations with the shock to the observed
        // variables, and that cross covariance is what a FAVAR is estimated to
        // measure.
        v = transition_residuals(state, a_mat, ns, p);
        v_sigma_inv = wishart(v, v_prior_scale, v_post_df);

        // Block 5: Draw the state transition ----
        //
        // The ns equations share their regressors, so the SUR design is
        // kron(X_a', I_ns) and its posterior precision collapses to
        // kron(X_a X_a', Q^-1). The identity holds for a full Q exactly as it
        // does for the diagonal one a dynamic factor model has -- it is a
        // statement about the regressors, not about the covariance.
        if (use_a)
        {
            fill_lagged_factors(x_a, state, ns, p);

            const arma::mat post_v =
                a_prior_vinv + arma::kron(x_a * arma::trans(x_a), v_sigma_inv);
            a = draw_normal_precision(
                post_v, a_prior_vinv * a_prior_mu +
                            arma::vectorise(v_sigma_inv * (state * arma::trans(x_a))));
            a_mat = arma::reshape(a, ns, ns * p);
        }

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            out.lambda.col(draw_pos) = arma::vectorise(lambda);
            out.factors.col(draw_pos) = arma::vectorise(factors);
            out.u_sigma_inv.col(draw_pos) = u_sigma_inv;
            out.v_sigma_inv.col(draw_pos) = arma::vectorise(v_sigma_inv);
            if (use_a)
            {
                out.a.col(draw_pos) = a;
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws FavarNormalWishartSampler::forecast(const FavarNormalWishartInput &input,
                                                  const FavarNormalWishartDraws &coefficients,
                                                  Reporter &reporter) const
{
    const int k = input.spec.k;
    const int n = input.spec.n_factors;
    const int n_obs = input.spec.n_obs_factors;
    const int ns = input.spec.n_state();
    const int p = input.spec.p;
    const int h = input.spec.h;
    const bool use_a = input.use_a();

    if (k <= 0 || n <= 0 || n_obs <= 0)
    {
        throw std::invalid_argument("a factor augmented VAR must have at least one observed "
                                    "series (k), one factor (n_factors) and one observed factor "
                                    "(n_obs_factors)");
    }
    if (h <= 0)
    {
        throw std::invalid_argument("forecast horizon (h) must be positive");
    }
    if (!coefficients.has_factors())
    {
        throw std::invalid_argument("posterior draws of the factors are missing; a factor "
                                    "augmented VAR forecasts by running the transition on from "
                                    "them");
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
        throw std::invalid_argument("the state has a transition of order " + std::to_string(p) +
                                    " but posterior draws of a are missing");
    }

    const int tt = static_cast<int>(input.train.periods(k));
    if (tt < p)
    {
        throw std::invalid_argument("a transition of order " + std::to_string(p) +
                                    " needs that many states to start from, and the sample has " +
                                    std::to_string(tt));
    }
    if (coefficients.factors.n_rows != static_cast<arma::uword>(n * tt))
    {
        throw std::invalid_argument(
            "posterior draws of the factors must have " + std::to_string(n * tt) + " rows for " +
            std::to_string(n) + " factors over " + std::to_string(tt) + " periods, got " +
            std::to_string(coefficients.factors.n_rows));
    }

    // The observed half of the history the transition starts from. A dynamic
    // factor model takes its whole starting state out of the posterior; here
    // half of it is data, and this is where that half comes from.
    const arma::mat obs = core::obs_factors_by_period(input.train);
    if (obs.n_cols != static_cast<arma::uword>(tt) ||
        obs.n_rows != static_cast<arma::uword>(n_obs))
    {
        throw std::invalid_argument("the observed factors must cover the whole sample, " +
                                    std::to_string(n_obs) + " of them over " +
                                    std::to_string(tt) + " periods");
    }

    const arma::uword draws = coefficients.iterations();
    const int width = k + n_obs;
    arma::mat fcst(h * width, draws);

    // The path carries the p states the transition needs before the first
    // horizon, so column p + i is horizon i and the lag lookup is one expression
    // at every horizon rather than a split between history and forecast.
    arma::mat path(ns, p + h);

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(draws));

        const arma::mat lambda = arma::reshape(coefficients.lambda.col(draw), k, ns);
        const arma::vec u_sd = 1.0 / arma::sqrt(coefficients.u_sigma_inv.col(draw));
        const arma::mat v_root =
            covariance_root(arma::reshape(coefficients.v_sigma_inv.col(draw), ns, ns));

        if (p > 0)
        {
            const arma::mat drawn = arma::reshape(coefficients.factors.col(draw), n, tt);
            path.submat(0, 0, n - 1, p - 1) = drawn.tail_cols(p);
            path.submat(n, 0, ns - 1, p - 1) = obs.tail_cols(p);
        }

        const arma::mat a_mat =
            use_a ? arma::reshape(coefficients.a.col(draw), ns, ns * p) : arma::mat();

        for (int i = 0; i < h; i++)
        {
            arma::vec s = v_root * arma::randn<arma::vec>(ns);
            for (int j = 1; j <= p; j++)
            {
                s += a_mat.cols((j - 1) * ns, j * ns - 1) * path.col(p + i - j);
            }
            path.col(p + i) = s;

            // The panel first, then the observed factors of the same horizon.
            // They are part of the state and so are forecast rather than
            // measured, which is why they are not read off the loadings.
            fcst.submat(i * width, draw, i * width + k - 1, draw) =
                lambda * s + u_sd % arma::randn<arma::vec>(k);
            fcst.submat(i * width + k, draw, (i + 1) * width - 1, draw) = s.tail(n_obs);
        }
    }

    reporter.finish();
    return ForecastDraws{fcst};
}

arma::mat
FavarNormalWishartSampler::log_likelihood(const FavarNormalWishartInput &input,
                                          const FavarNormalWishartDraws &coefficients) const
{
    const int k = input.spec.k;
    const int n = input.spec.n_factors;
    const int ns = input.spec.n_state();

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
        throw std::invalid_argument("posterior draws of the factors are missing; the likelihood "
                                    "of a factor augmented VAR is conditional on them");
    }
    if (coefficients.lambda.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of lambda are missing");
    }

    const int tt = static_cast<int>(input.train.periods(k));
    const arma::mat x_t = response_by_period(input.train, k, tt);
    const arma::mat obs = core::obs_factors_by_period(input.train);
    const arma::uword draws = coefficients.iterations();

    arma::mat loglik(draws, tt);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        const arma::mat lambda = arma::reshape(coefficients.lambda.col(draw), k, ns);
        const arma::mat factors = arma::reshape(coefficients.factors.col(draw), n, tt);
        const arma::vec u_sigma_inv = coefficients.u_sigma_inv.col(draw);

        // The panel is scored against the whole state, drawn half and observed
        // half alike. The observed factors are data the state is conditioned on
        // rather than a further n_obs columns of fit, so they contribute to what
        // the panel is explained by and not to the count of what is explained.
        const arma::mat state = stacked_state(factors, obs);

        // R is diagonal, so the determinant term is a sum of logs and the
        // quadratic form is a weighted sum of squares -- no k x k anything.
        const double part_b = arma::accu(arma::log(u_sigma_inv)) / 2;
        const arma::mat u = x_t - lambda * state;

        for (int i = 0; i < tt; i++)
        {
            const double part_c = -arma::dot(u_sigma_inv, arma::square(u.col(i))) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
