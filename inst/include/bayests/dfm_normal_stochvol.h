// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_DFM_NORMAL_STOCHVOL_H
#define BAYESTS_DFM_NORMAL_STOCHVOL_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// Dynamic factor model with constant loadings, a constant factor transition and
/// stochastic volatility in both error terms.
///
///     x_t = Lambda f_t + u_t,                  u_t ~ N(0, U_t),  U_t diagonal,
///     f_t = sum_{j=1..p} A_j f_{t-j} + v_t,    v_t ~ N(0, V_t),  V_t diagonal,
///
/// with U_t = diag(exp(h^u_t)), V_t = diag(exp(h^v_t)) and every element of both
/// log-volatilities a random walk of its own. DfmNormalGamma with the two gamma
/// priors on the precisions replaced by two stochastic volatility blocks, and
/// nothing else changed: `Normal` still says the loadings and the transition are
/// constant, exactly as it does in VarNormalStochvol against VarTvpStochvol.
///
/// Seven Gibbs blocks -- the factor path, the loadings, the two log-volatilities
/// with their state equations, and the transition -- against DfmNormalGamma's
/// five.
///
/// What the extra volatility changes, block by block, is more than the shape of
/// the output:
///
///   - The factor path is still one Gaussian vector of block banded precision and
///     is still drawn by `chan_jeliazkov_2009`, which already takes a covariance
///     per period in both equations. Nothing new was needed there. It does need
///     the period indexing of the transition covariances to be right, and that is
///     off by one from this model's; `draw_factor_path_sv` is where that is dealt
///     with and says why.
///   - The loadings are a weighted regression rather than a scaled one. Row i of
///     Lambda regresses on the factors with weights exp(-h^u_{it}), so the
///     periods in which series i was quiet identify it and the periods in which
///     it was wild mostly do not. That reweighting, not the wider output, is the
///     reason to prefer this model to DfmNormalGamma on a sample that spans a
///     change in volatility.
///   - The transition loses the Kronecker collapse. DfmNormalGamma's posterior
///     precision is kron(X X', V^-1) because V is the same in every period;
///     sum_t kron(x_t x_t', V_t^-1) does not factor. What survives is that V_t is
///     diagonal, so the N equations of the transition are conditionally
///     independent and each contributes its own weighted cross-product. See the
///     source: no (tt N) x (N^2 p) matrix is built here either.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG -- the same contract DfmNormalGamma keeps, and what lets one
/// object serve the command line and an embedded caller such as an R package.
///
/// Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). Bayesian
/// econometric methods (2nd ed.). Cambridge: Cambridge University Press.
///
/// Omori, Y., Chib, S., Shephard, N., & Nakajima, J. (2007). Stochastic
/// volatility with leverage: Fast and efficient likelihood inference. Journal of
/// Econometrics, 140(2), 425-449.
class DfmNormalStochvolSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    DfmNormalStochvolDraws draw_coefficients(const DfmNormalStochvolInput &input,
                                             Reporter &reporter) const;

    /// Simulates one forecast path of the observed series per posterior draw.
    ///
    /// As in DfmNormalGamma, no out-of-sample regressor matrix is needed: the
    /// path is the transition run forward from the last p drawn factors, with an
    /// innovation drawn at each step, and the observed series read off the
    /// loadings. `spec.h` is the horizon.
    ///
    /// Both volatilities are held at their terminal value over the horizon rather
    /// than simulated forward. That is the convention every stochastic volatility
    /// model here follows -- VarNormalStochvol reads its precision at period
    /// tt - 1 and uses it at every horizon -- and it is what the posterior
    /// supports: the variance of the log-volatility innovations is a state of the
    /// chain, not something the draws carry, so there is nothing to extrapolate
    /// the random walk with. The predictive intervals are therefore conditional on
    /// the volatility that prevailed at the end of the sample, which understates
    /// them by however much the volatility might still move.
    ///
    /// Requires `draws.factors`, both precisions, and `draws.a` when the
    /// transition has an order. `u_sigma_inv` and `v_sigma_inv` may hold either
    /// the whole path or the terminal period alone, which is what lets a host
    /// read one period out of a file instead of all of it.
    ///
    /// Throws std::invalid_argument if a required draw is missing or if spec.h is
    /// zero.
    ForecastDraws forecast(const DfmNormalStochvolInput &input,
                           const DfmNormalStochvolDraws &draws, Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior draw,
    /// one column per observation, as expected by WAIC and PSIS-LOO.
    ///
    /// Conditional on the drawn factor path and on the drawn volatility: what is
    /// evaluated is p(x_t | f_t, Lambda, U_t). The reasoning DfmNormalGamma gives
    /// for conditioning on the factors applies to the volatility as well -- both
    /// are states this sampler draws and stores, and integrating either out would
    /// cost a filtering pass per draw and answer a different question.
    ///
    /// Requires `draws.factors` and the whole path of `u_sigma_inv`.
    arma::mat log_likelihood(const DfmNormalStochvolInput &input,
                             const DfmNormalStochvolDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_DFM_NORMAL_STOCHVOL_H
