// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_DFM_TVP_STOCHVOL_H
#define BAYESTS_DFM_TVP_STOCHVOL_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// Dynamic factor model whose loadings and factor transition follow random walks
/// and whose two error terms carry stochastic volatility.
///
///     x_t = Lambda_t f_t + u_t,                  u_t ~ N(0, U_t),
///     f_t = sum_{j=1..p} A_{j,t} f_{t-j} + v_t,  v_t ~ N(0, V_t),
///
/// with U_t = diag(exp(h^u_t)) and V_t = diag(exp(h^v_t)). The widest model here:
/// DfmTvpGamma's coefficients over DfmNormalStochvol's errors, and nothing in it
/// is held fixed but the normalisation. Nine Gibbs blocks -- the factor path, the
/// loading path with its state variance and initial state, the two
/// log-volatilities with theirs, and the transition path with the same pair.
///
/// Why both halves. A series' exposure to the common component is not a constant
/// of nature, and neither is the scale of the shocks it is measured against; a
/// model with only one of the two has to explain the other with what it has. A
/// series whose loading fell looks like a series whose idiosyncratic variance
/// rose, and a period of common turbulence looks like a transition that changed.
/// Carrying both is what lets the data say which, and it is what Del Negro and
/// Otrok (2008) put together for the same reason.
///
/// The identifying block still does not drift. Lambda's leading N x N block stays
/// unit lower triangular in every period, exactly as in DfmTvpGamma, because only
/// the product Lambda_t f_t is identified.
///
/// What this model costs, and where. The factor path needs no new algorithm --
/// `chan_jeliazkov_2009` takes a measurement matrix, a measurement covariance, a
/// transition and a transition covariance per period, and this is the one model
/// that stacks all four -- but a per-period measurement means its assembly forms
/// Z' U^-1 Z once per period rather than once, which with many observed series is
/// the dominant cost. `draw_factor_path()` carries the shift conventions that
/// arrangement needs, all of them in one place.
///
/// The loading paths are drawn row by row, each against its own series'
/// volatility: row i is weighted by exp(-h^u_t(i)) period by period, so the
/// periods in which series i was quiet identify its loadings and the periods in
/// which it was wild largely do not. That interaction is the one thing here that
/// neither DfmTvpGamma nor DfmNormalStochvol has on its own.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
///
/// Del Negro, M., & Otrok, C. (2008). Dynamic factor models with time-varying
/// parameters: measuring changes in international business cycles. Federal
/// Reserve Bank of New York Staff Report No. 326.
///
/// Omori, Y., Chib, S., Shephard, N., & Nakajima, J. (2007). Stochastic
/// volatility with leverage. Fast and efficient likelihood inference. Journal of
/// Econometrics, 140(2), 425-449.
class DfmTvpStochvolSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    DfmTvpStochvolDraws draw_coefficients(const DfmTvpStochvolInput &input,
                                          Reporter &reporter) const;

    /// Simulates one forecast path of the observed series per posterior draw,
    /// holding the loadings, the transition and both volatilities at their last
    /// in-sample values.
    ///
    /// `draws.lambda` and `draws.a` are expected to carry that period alone, one
    /// column per draw; the two precisions may carry either the whole path or the
    /// terminal block, since the last block of a path is what is read either way.
    ///
    /// As in every dynamic factor model here, no out-of-sample regressor matrix
    /// is needed: the path is the transition run forward from the last p drawn
    /// factors, with an innovation drawn at each step.
    ///
    /// Requires `draws.factors`, and `draws.a` when the transition has an order.
    ForecastDraws forecast(const DfmTvpStochvolInput &input, const DfmTvpStochvolDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior draw,
    /// one column per observation, as expected by WAIC and PSIS-LOO.
    ///
    /// Every period is scored under its own loadings *and* its own precision, so
    /// both paths are wanted whole here; the determinant term moves inside the
    /// period loop for the second reason and the fitted value for the first.
    /// Conditional on the drawn factor path, for the reason DfmNormalGamma's
    /// log_likelihood() sets out at length.
    ///
    /// Requires `draws.factors`.
    arma::mat log_likelihood(const DfmTvpStochvolInput &input,
                             const DfmTvpStochvolDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_DFM_TVP_STOCHVOL_H
