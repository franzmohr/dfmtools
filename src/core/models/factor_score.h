// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_FACTOR_SCORE_H
#define BAYESTS_CORE_MODELS_FACTOR_SCORE_H

#include "bayests/arma.h"
#include "bayests/spec.h"

#include "core/models/predictive_score.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace bayests::core
{

/// Scoring a factor model, which is a different problem from scoring a VAR.
///
/// Every other model here reaches its realised history through its regressors:
/// the lag blocks of a scored period carry what was realised, so the regressors
/// do not depend on the draw, and the score is the model's own pointwise log
/// likelihood on a different sample. A factor model's history reaches the
/// density through the **factors**, which are latent. Its in-sample likelihood
/// conditions on the factors the sampler drew, and no such draw exists for a
/// period outside the sample.
///
/// Scoring one as the columns of `/posterior/forecast/loglik` are defined --
/// each conditioning on the observations realised before it -- therefore means
/// filtering. At every scored period the realised `y` updates the distribution
/// of the factors before the next period is predicted, and the column is that
/// period's prediction error decomposition:
///
///     x_{i|i-1} = T x_{i-1|i-1},        P_{i|i-1} = T P_{i-1|i-1} T' + Q
///     v_i       = y*_i - Z x_{i|i-1},   F_i       = Z P_{i|i-1} Z' + R
///     log p     = -k/2 log 2pi - log|F_i| / 2 - v_i' F_i^-1 v_i / 2
///     x_{i|i}   = x_{i|i-1} + K v_i,    P_{i|i}   = P_{i|i-1} - K Z P_{i|i-1}
///
/// with K = P_{i|i-1} Z' F_i^-1. Summing those columns over the periods is the
/// joint log density of the realised stretch, factorised, exactly as it is for a
/// VAR -- which is the whole reason for filtering rather than simulating the
/// factors forward. Simulating them, as the forecast does, would give the
/// marginal density of each horizon instead: a defensible quantity, and not the
/// one this dataset holds for every other algorithm.
///
/// **The filter starts from the drawn factors, with no uncertainty.** Conditional
/// on a draw, the factors at the end of the sample are known -- they are part of
/// what was drawn -- so P at the start is zero, and averaging exp() over draws
/// integrates them out with everything else the posterior carries. That is the
/// same convention the forecast uses when it runs the transition on from the
/// drawn path.

/// What a factor model is at one scored period. A model whose loadings or
/// variances stand still fills the same values at every period; one whose drift
/// fills what that period's states give.
struct FactorPeriod
{
    /// The factors the transition needs before the first scored period, `n` rows
    /// by `p` columns with the most recent last. Read at `i == 0` only.
    arma::mat start;

    arma::mat lambda;     ///< k x n, the loadings.
    arma::mat transition; ///< n x (n*p), [A_1 .. A_p]; empty where the factors have no dynamics.
    arma::vec v_var;      ///< n, the variances of the factor innovations.
    arma::vec u_var;      ///< k, the idiosyncratic variances.
};

/// The score of a factor model, draws x scored periods.
///
/// `step(draw, i, out)` is called once per scored period, i = 0 included, and
/// fills `out` with that period's parameters -- and, at i == 0, with the factors
/// the filter starts from. It draws nothing: a filter is a recursion, not a
/// simulation, so unlike a forecast the score of a factor model is computed
/// rather than drawn and repeats without a seed.
template <typename Step>
inline arma::mat score_factor_forecast(const VarSpec &spec, const arma::mat &realised,
                                       const arma::uword draws, Step &&step)
{
    require_scorable(spec, "a factor model");

    const arma::uword periods = scored_horizons(realised, spec);
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword n = static_cast<arma::uword>(spec.n_factors);
    const arma::uword p = spec.p > 0 ? static_cast<arma::uword>(spec.p) : 0;

    if (k == 0 || n == 0)
    {
        throw std::invalid_argument(
            "scoring a factor model needs at least one observed series (k) and one factor");
    }
    if (realised.n_cols != k)
    {
        throw std::invalid_argument(
            "the realised values must have k = " + std::to_string(spec.k) + " columns, got " +
            std::to_string(realised.n_cols));
    }

    // The state stacks the p factors the transition reaches back over. With no
    // dynamics there is one block and the transition is zero: the factors are
    // then white noise and the filter still says the right thing, the prediction
    // being the mean of that noise.
    const arma::uword blocks = p > 0 ? p : 1;
    const arma::uword n_state = n * blocks;

    arma::mat loglik(draws, periods);
    const double constant = -static_cast<double>(k) * std::log(2 * arma::datum::pi) / 2;

    FactorPeriod current;
    arma::mat transition(n_state, n_state, arma::fill::zeros);
    arma::mat state_noise(n_state, n_state, arma::fill::zeros);
    arma::mat observation(k, n_state, arma::fill::zeros);
    arma::vec state(n_state);
    arma::mat variance(n_state, n_state);

    // The shift that carries yesterday's blocks down one, which never changes.
    if (blocks > 1)
    {
        transition.submat(n, 0, n_state - 1, n_state - n - 1) =
            arma::eye<arma::mat>(n_state - n, n_state - n);
    }

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        for (arma::uword i = 0; i < periods; i++)
        {
            step(draw, static_cast<int>(i), current);

            if (i == 0)
            {
                if (p > 0)
                {
                    if (current.start.n_rows != n || current.start.n_cols != p)
                    {
                        throw std::invalid_argument(
                            "the filter starts from the last " + std::to_string(p) +
                            " periods of the drawn factors, " + std::to_string(n) +
                            " per period");
                    }
                    // Most recent first, which is the order the state stacks in.
                    for (arma::uword j = 0; j < p; j++)
                    {
                        state.subvec(j * n, (j + 1) * n - 1) = current.start.col(p - 1 - j);
                    }
                }
                else
                {
                    state.zeros();
                }

                // Conditional on the draw the factors are known, so the filter
                // starts certain of them.
                variance.zeros();
            }

            if (p > 0)
            {
                transition.submat(0, 0, n - 1, n_state - 1) = current.transition;
            }
            state_noise.submat(0, 0, n - 1, n - 1) = arma::diagmat(current.v_var);
            observation.submat(0, 0, k - 1, n - 1) = current.lambda;

            // Predict.
            state = transition * state;
            variance = transition * variance * arma::trans(transition) + state_noise;

            // The one step ahead forecast of the observation, and its error.
            const arma::vec innovation = arma::trans(realised.row(i)) - observation * state;
            const arma::mat gain = variance * arma::trans(observation);
            arma::mat scale = observation * gain + arma::diagmat(current.u_var);
            scale = (scale + arma::trans(scale)) / 2;

            // Factorised rather than put through arma::log_det(), which is an LU
            // and decides definiteness from a pivot count. F is symmetric
            // positive definite by construction, so a Cholesky is the test that
            // matches it -- it fails exactly when the matrix is not one -- and it
            // carries the determinant and both solves below with it.
            arma::mat chol_scale;
            if (!arma::chol(chol_scale, scale, "lower"))
            {
                throw std::runtime_error(
                    "the one step ahead forecast variance of a scored period is not positive "
                    "definite, which no set of variances this filter was given can make it");
            }
            const double log_det = 2.0 * arma::accu(arma::log(chol_scale.diag()));

            const arma::vec weighted =
                arma::solve(arma::trimatu(arma::trans(chol_scale)),
                            arma::solve(arma::trimatl(chol_scale), innovation));
            loglik(draw, i) =
                constant - log_det / 2 - arma::dot(innovation, weighted) / 2;

            // Update, so that the next period conditions on what this one
            // realised, in the **Joseph form**: P = (I - K Z) P (I - K Z)' +
            // K R K', a sum of two positive semi-definite terms and so one
            // whatever rounding does to it.
            //
            // The short form P - K F K' the prediction error decomposition is
            // written against is a *difference* of two of them, and it is not
            // safe here for a reason particular to this filter: `state_noise`
            // carries diag(v_var) in the leading n x n block and zeros
            // elsewhere, so Q is singular for any transition of order above one,
            // and the lagged blocks of P are never refreshed by it -- they only
            // shift down through T. Error there accumulates instead of being
            // flooded out, P drifts indefinite, and the Cholesky above then
            // rejects an F that is mathematically fine. That took a release to
            // find, because whether the drift crosses zero depends on which BLAS
            // rounds which way: it passed on two toolchains and failed on a
            // third.
            const arma::mat kalman =
                arma::trans(arma::solve(arma::trimatu(arma::trans(chol_scale)),
                                        arma::solve(arma::trimatl(chol_scale),
                                                    arma::trans(gain))));
            state += kalman * innovation;

            const arma::mat spread =
                arma::eye<arma::mat>(n_state, n_state) - kalman * observation;
            variance = spread * variance * arma::trans(spread) +
                       kalman * arma::diagmat(current.u_var) * arma::trans(kalman);
            variance = (variance + arma::trans(variance)) / 2;
        }
    }

    return loglik;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_FACTOR_SCORE_H
