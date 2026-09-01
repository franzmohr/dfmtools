// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_DFM_NORMAL_GAMMA_H
#define BAYESTS_DFM_NORMAL_GAMMA_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// Dynamic factor model with a normal prior on the loadings and on the factor
/// transition, and independent gamma priors on both error precisions.
///
///     x_t = Lambda f_t + u_t,                  u_t ~ N(0, U),  U diagonal,
///     f_t = sum_{j=1..p} A_j f_{t-j} + v_t,    v_t ~ N(0, V),  V diagonal,
///
/// for M observed series and N unobserved factors, after Chan, Koop, Poirier and
/// Tobias (2019). Five Gibbs blocks: the factor path, the loadings, the two
/// precisions and the transition.
///
/// What makes this the first model here that is not a regression. The factors
/// are unobserved, so the sampler draws a whole tt-period path of them every
/// iteration and carries it in the posterior alongside the parameters; and the
/// loadings are identified only up to a rotation and a scale, so the leading
/// N x N block of Lambda is fixed unit lower triangular rather than drawn. That
/// second point is why the loadings are drawn equation by equation instead of as
/// one vector: row i of Lambda has min(i, N) free elements against different
/// regressors, so the k rows do not share a design matrix the way the equations
/// of a VAR do.
///
/// The factor path is drawn whole, from its precision, by
/// `chan_jeliazkov_2009`. That is exactly the object it was written for: the
/// transition is order p, so the path's posterior precision is block banded of
/// bandwidth p, and the alternative -- forming the (tt N) square matrix the
/// reference implementation builds and factorising it -- is O(tt^3 N^3) against
/// O(tt N^3). Values of the factors before the sample are zero rather than
/// drawn, which is what pins down the first p transitions; the prior on the
/// first p states that this implies is derived in the source.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
///
/// Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). Bayesian
/// econometric methods (2nd ed.). Cambridge: Cambridge University Press.
class DfmNormalGammaSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    DfmNormalGammaDraws draw_coefficients(const DfmNormalGammaInput &input,
                                          Reporter &reporter) const;

    /// Simulates one forecast path of the observed series per posterior draw.
    ///
    /// Unlike every other forecast here, this one needs no out-of-sample
    /// regressor matrix. A DFM has no regressors: the path is the transition run
    /// forward from the last p drawn factors, with an innovation drawn at each
    /// step, and the observed series read off the loadings. `spec.h` is the only
    /// thing that says how far.
    ///
    /// Requires `draws.factors`, and `draws.a` when the transition has an order.
    /// Throws std::invalid_argument if either is missing or if spec.h is zero.
    ForecastDraws forecast(const DfmNormalGammaInput &input, const DfmNormalGammaDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior draw,
    /// one column per observation, as expected by WAIC and PSIS-LOO.
    ///
    /// Conditional on the drawn factor path: what is evaluated is
    /// p(x_t | f_t, Lambda, U), the measurement density at the factors that were
    /// stored beside the parameters. That is the conditional -- not the marginal
    /// -- likelihood, and the distinction matters for what the resulting
    /// information criterion means: conditional WAIC assesses the fit of a model
    /// that treats the states as parameters, which is the quantity a sampler that
    /// draws them can report without a filtering pass per draw. Integrating the
    /// factors out would be a different number and would cost a Kalman filter
    /// over the sample for every draw.
    ///
    /// Requires `draws.factors`.
    arma::mat log_likelihood(const DfmNormalGammaInput &input,
                             const DfmNormalGammaDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_DFM_NORMAL_GAMMA_H
