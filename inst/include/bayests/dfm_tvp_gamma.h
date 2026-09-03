// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_DFM_TVP_GAMMA_H
#define BAYESTS_DFM_TVP_GAMMA_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// Dynamic factor model whose loadings and factor transition follow random
/// walks, with independent gamma priors on both error precisions.
///
///     x_t = Lambda_t f_t + u_t,                  u_t ~ N(0, U),  U diagonal,
///     f_t = sum_{j=1..p} A_{j,t} f_{t-j} + v_t,  v_t ~ N(0, V),  V diagonal,
///
/// DfmNormalGamma with its two coefficient blocks turned into state paths and
/// nothing else changed: `Tvp` names the coefficients, as it does in VarTvpGamma
/// against VarNormalGamma, and the errors stay homoskedastic. Seven Gibbs blocks
/// against five -- the factor path, the loading path with its state variance and
/// its initial state, the two precisions, and the transition path with the same
/// pair.
///
/// What the drift is for. A loading is a series' exposure to the common factor,
/// and the assumption that it held over the whole sample is the one a factor
/// model makes most often and defends least: a series can enter or leave the
/// common component -- a sector reorganised, a country's trade opening -- without
/// anything about the factor itself changing. A constant-loading model has
/// nowhere to put that except the idiosyncratic variance, which then absorbs it
/// as noise the series is credited with throughout, including in the periods
/// where the exposure did hold. Drift in the transition is the other half: the
/// persistence of the common component is what a forecast from it runs on, and
/// it is not a constant of nature either.
///
/// The identifying block does not drift. Lambda's leading N x N block stays unit
/// lower triangular in every period, exactly as in DfmNormalGamma, because only
/// the product Lambda_t f_t is identified: letting the block move would let the
/// rotation and the scale of the factors wander over the sample, and a loading
/// path would then describe the normalisation as much as the exposure it is read
/// as. Every free element drifts; none of the fixed ones does.
///
/// Three things in the numerics are worth knowing about.
///
///   - The factor path needed no new algorithm. `chan_jeliazkov_2009` already
///     took a measurement matrix and a transition per period, so a drifting
///     Lambda and a drifting A reach it as stacks. What is lost is the shortcut
///     it takes when the measurement is the same in every period: a dynamic
///     factor model is worth having when M is large, and the assembly's
///     Z'U^-1 Z is then the dominant cost, paid tt times here instead of once.
///     That is the price of the model, not an implementation choice.
///   - The loading paths are drawn row by row, which is not an economy but the
///     shape of the problem: row i has min(i, N) free elements against the first
///     min(i, N) factors, so the rows share no design, and given the factors and
///     a diagonal U they are conditionally independent. Each row is one call to
///     `kalman_durbin_koopman_2002` over a state of its own width.
///   - The transition path is drawn as one state of N^2 p elements, from the SUR
///     design kron(x_t', I_N) built per period. The Kronecker identity
///     DfmNormalGamma leans on -- the N equations sharing their regressors, so
///     the posterior precision collapses to kron(X X', V^-1) -- is a statement
///     about a single coefficient vector and does not survive the coefficients
///     becoming a path.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
///
/// Del Negro, M., & Otrok, C. (2008). Dynamic factor models with time-varying
/// parameters: measuring changes in international business cycles. Federal
/// Reserve Bank of New York Staff Report No. 326.
class DfmTvpGammaSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    DfmTvpGammaDraws draw_coefficients(const DfmTvpGammaInput &input, Reporter &reporter) const;

    /// Simulates one forecast path of the observed series per posterior draw,
    /// holding the loadings and the transition at their last in-sample values:
    /// `draws.lambda` and `draws.a` are expected to carry that period alone, one
    /// column per draw.
    ///
    /// As in DfmNormalGamma, this needs no out-of-sample regressor matrix. A DFM
    /// has no regressors: the path is the transition run forward from the last p
    /// drawn factors, with an innovation drawn at each step, and the observed
    /// series read off the loadings.
    ///
    /// Requires `draws.factors`, and `draws.a` when the transition has an order.
    /// Throws std::invalid_argument if either is missing or if spec.h is zero.
    ForecastDraws forecast(const DfmTvpGammaInput &input, const DfmTvpGammaDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior draw,
    /// one column per observation, as expected by WAIC and PSIS-LOO.
    ///
    /// `draws.lambda` carries the whole loading path: every period is scored
    /// under its own Lambda_t, against the single idiosyncratic precision this
    /// model has. Conditional on the drawn factor path, for the reason
    /// DfmNormalGamma's log_likelihood() sets out at length -- integrating the
    /// factors out would be a different number and would cost a filtering pass
    /// per draw.
    ///
    /// Requires `draws.factors`.
    arma::mat log_likelihood(const DfmTvpGammaInput &input,
                             const DfmTvpGammaDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_DFM_TVP_GAMMA_H
