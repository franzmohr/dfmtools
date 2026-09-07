// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_FAVAR_NORMAL_WISHART_H
#define BAYESTS_FAVAR_NORMAL_WISHART_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// Factor augmented VAR with a normal prior on the loadings and on the state
/// transition, independent gamma priors on the idiosyncratic precisions and a
/// Wishart prior on the precision of the state innovations.
///
///     x_t = Lambda_f f_t + Lambda_y y_t + e_t,   e_t ~ N(0, R),  R diagonal,
///     s_t = sum_{j=1..p} Phi_j s_{t-j} + v_t,    v_t ~ N(0, Q),
///
/// with s_t = (f_t', y_t')', for M observed series in the panel, N unobserved
/// factors and Ky observed ones, after Bernanke, Boivin and Eliasz (2005). Five
/// Gibbs blocks: the factor path, the loadings, the idiosyncratic precisions,
/// the state innovation precision and the transition.
///
/// What separates this from the dynamic factor models beside it, and it is one
/// thing. A DFM's state is unobserved throughout, so it is drawn whole. Half of
/// this state is data -- the observed factors are the variables the model is
/// about, a policy rate or output, and the panel is there to measure the common
/// component they move with -- so the path draw *conditions* on that half rather
/// than drawing it. They sit in the state vector rather than beside it because
/// the transition is a VAR over both blocks jointly: the factors respond to the
/// policy variable and the policy variable responds to the factors, and that
/// coupling is the model.
///
/// Conditioning, not approximating. Adding the observed factors to the
/// measurement equation with a small error variance is the usual shortcut, and
/// it is not available to a precision based sampler at all -- observing
/// something exactly is infinite precision, which the band cannot hold. What
/// `chan_jeliazkov_2009_conditional` does instead is exact: the assembled
/// precision is partitioned into the drawn rows and the observed ones, and
/// K_FF f = b_F - K_FY y is the same band one block size narrower. Nothing is
/// dropped and nothing is approximated, and the information the observed
/// factors' own equations carry about the lagged factors -- which a sampler that
/// drew the factor block alone would throw away -- is in there through the cross
/// terms.
///
/// Why the error precision is a Wishart here and a gamma in every DFM. R must
/// stay diagonal: a factor model whose idiosyncratic errors may correlate has
/// nothing left for the factors to explain, which is why no `Dfm*` offers the
/// choice. Q is a different object. It is the innovation covariance of a VAR;
/// its observed block is an ordinary VAR covariance and its cross block is the
/// correlation between the factor innovations and the shock to the observed
/// variables, which is the one thing a FAVAR exists to measure. So the third
/// part of the name refers to Q, and R is gamma-diagonal in every member of this
/// family -- the one place the naming rule differs from the DFM row above it.
///
/// The identification is the DFM's restriction made stronger, and it has to be.
/// Lambda's leading N x N block is the *identity* in the factor columns -- not
/// the unit lower triangle a DFM uses -- and zero in the observed ones, so the
/// first N series of the panel are the factors plus idiosyncratic noise and
/// carry no free loading at all.
///
/// The stronger form is forced by the Wishart. A rotation F -> C F leaves the
/// measurement unchanged if Lambda_f absorbs it, and a DFM rules C out with two
/// restrictions working together: a unit lower triangular block and a diagonal
/// V, which jointly admit only C = I by the uniqueness of an LDL factorisation.
/// A FAVAR has no diagonal V to offer -- Q is free, that being the model -- so a
/// unit lower triangular block on its own would leave every unit lower
/// triangular C admissible and the loadings free to wander along a ridge. The
/// identity block rules out C by itself, and the zero block rules out the other
/// transformation the state admits, F -> F + D Y. Together they pin everything,
/// which is why Bernanke, Boivin and Eliasz identify the model this way.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
///
/// Bernanke, B. S., Boivin, J., & Eliasz, P. (2005). Measuring the effects of
/// monetary policy: a factor-augmented vector autoregressive (FAVAR) approach.
/// Quarterly Journal of Economics, 120(1), 387-422.
class FavarNormalWishartSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    FavarNormalWishartDraws draw_coefficients(const FavarNormalWishartInput &input,
                                              Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw.
    ///
    /// As in a DFM this needs no out-of-sample regressor matrix: the path is the
    /// transition run forward from the last p states, with an innovation drawn
    /// at each step, and the panel read off the loadings. Unlike a DFM the last
    /// p states are half drawn and half data -- the observed factors come from
    /// `input.train.f_obs`, which is where the model's own history of them is.
    ///
    /// The result is the one forecast here that is wider than `k`: each horizon
    /// carries the k panel series followed by the n_obs_factors observed ones,
    /// which are what a FAVAR is forecast for and have no other dataset to go
    /// in. See ForecastDraws.
    ///
    /// Requires `draws.factors`, and `draws.a` when the state has a transition.
    /// Throws std::invalid_argument if either is missing or if spec.h is zero.
    ForecastDraws forecast(const FavarNormalWishartInput &input,
                           const FavarNormalWishartDraws &draws, Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior draw,
    /// one column per observation, as expected by WAIC and PSIS-LOO.
    ///
    /// The measurement density of the panel, p(x_t | s_t, Lambda, R),
    /// conditional on the drawn factor path and on the observed factors. That is
    /// the conditional -- not the marginal -- likelihood, for the reason
    /// DfmNormalGamma::log_likelihood() sets out at length, and it scores the
    /// panel alone: the observed factors are data the state is conditioned on
    /// rather than a further k columns of fit.
    ///
    /// Requires `draws.factors`.
    arma::mat log_likelihood(const FavarNormalWishartInput &input,
                             const FavarNormalWishartDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_FAVAR_NORMAL_WISHART_H
