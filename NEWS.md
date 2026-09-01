# dfmtools (development version)

* **Stochastic volatility in both error terms**, with
`create_dfmodel(error = "sv")`. The measurement equation becomes
$x_t = \lambda f_t + u_t$ with $u_t \sim N(0, U_t)$ and the transition
$f_t = \sum_i A_i f_{t-i} + v_t$ with $v_t \sim N(0, V_t)$, both covariances
diagonal and the log-volatility of every element a random walk of its own, drawn
with the mixture approximation of Omori et al. (2007). The loadings and the
transition stay constant, which is what the sampler's name says: the new
algorithm is `DfmNormalStochvol`, beside the existing `DfmNormalGamma`.

    Both placements are there because neither substitutes for the other.
    Volatility in $u_t$ reweights the series that identify the factors, so a
    series that was noisy early and quiet later stops contributing on the same
    terms throughout, and a single wild observation is absorbed where it happened
    rather than dragged into the factor. Volatility in $v_t$ is the common
    component's own, and it is what keeps the $M$ idiosyncratic variances from
    jointly absorbing a shock that every series felt at once — the factor is
    otherwise flattest exactly when it should move most.

    What changes for a caller of the existing model: nothing.
    `error = "gamma"` remains the default and its draws are unchanged. For
    `error = "sv"`, arguments `u` and `v` of `add_priors()` take the six-element
    stochastic volatility specification `bvartools` uses for its VAR and VEC
    models — `mu`, `v_i`, `shape`, `rate`, `state_variance`, `offset` — rather
    than `shape` and `rate` alone, and there are no defaults for it, so a call
    that forgets is told which elements are missing. `add_initial_values()`
    returns a log-volatility path per error term, `u_h` and `v_h`, in place of the
    precision matrices `uinv` and `vinv`. In the posterior, `u_sigma_inv` and
    `v_sigma_inv` widen from one number per series to a whole path — $M \times T$
    and $N \times T$ columns per draw — so the variance of series $i$ in period
    $t$ is `matrix(1 / draw, m, tt)[i, t]`. Forecasts hold both volatilities at
    their last in-sample value, as every stochastic volatility model in this
    family does.

    Note that `u` describes $M$ observed series and `v` describes $N$ factors, so
    the two specifications are the same shape but not the same length. Filling one
    from the other gives a well-formed list of the wrong width; the sampler
    rejects it by name.

* The vendored BayesTS core was refreshed to pick the new sampler up, and grew
five files: the sampler and its declaration, and the stochastic volatility
mixture it draws with. *Draws are unchanged* for the existing model — of the
shared files that moved, one gained a diagonal fast path that returns the same
numbers a dense Cholesky inverse of a diagonal matrix does, one gained a second
accepted argument shape with the old one bit-identical to before, and one moved
a block of checks between functions without altering them. BayesTS's own golden
fingerprint harness, 145 tests, passes.

# dfmtools 0.1.0

* Dynamic factor models now live here rather than in bvartools, and this package
is where they are maintained. The functions that moved -- `create_dfmodel`,
`add_priors.dfmodel`, `add_initial_values.dfmodel` and the `bem_dfmdata` data
set -- arrive in the state bvartools had them in, which was ahead of the copies
this package carried from 2024.

* **Posterior simulation is no longer implemented in R.** `dfmpost()` is
replaced by the vendored core layer of
[BayesTS](https://github.com/franzmohr/BayesTS), reached through
`add_posterior_coefficients()`. That is the arrangement bvartools uses for its
VAR and VEC models: one implementation of each sampler, upstream, with the R
package as a translation layer over it. `src/core/VENDORED.md` says which files
were copied, how the set is computed, and what was modified on the way in.

* **Draws change, and the previous ones were wrong outside a narrow case.**
`dfmpost()` was correct at one factor and a transition of order at most two.
Three defects put the limit there, and none is reproduced:

    * The loadings were drawn from the wrong distribution whenever there was more
      than one factor. `.post_lambda` drew with `solve(chol(K, "lower"), z)`,
      whose covariance is `(L'L)^-1` rather than the intended `K^-1`; the two
      agree only when the block is 1x1, which for a row of the loading matrix
      means one factor.
    * The transition's regressor matrix laid its lag blocks on top of one
      another for more than one factor. The block for lag `i` occupies rows
      `(i - 1) * n + 1:n`; `dfmpost()` wrote rows `(i - 1) + 1:n`.
    * `bvartools::generate_lower_block_diagonal()`, which built the transition's
      band, wrote past the end of the matrix for `p >= 3` and dropped a
      coefficient block from the last columns. The core builds the equivalent
      structure itself.

* Two smaller changes are deliberate rather than corrective. The free loadings
are ordered row by row wherever they appear as a vector -- the starting value
and both halves of the prior -- which is the order the equation-by-equation draw
consumes them in; `dfmpost()` stored them column-major but sliced the prior
precision row-major, a mismatch invisible only because `add_priors` builds that
precision as a scalar diagonal. And the posterior reports the whole `M x N`
loading matrix, the fixed ones and zeros of the identifying block included,
rather than the free elements alone, so a draw is reshaped rather than unpacked.

* The factor path is part of the posterior, in element `factors`, one whole
`T`-period path per draw. The factors are unobserved, so neither the forecast nor
the log-likelihood can be recomputed from the parameters and both read it back.

* Added `add_posterior_forecasts.dfmodel()`. A dynamic factor model needs no
out-of-sample regressors: each draw's path is the transition run forward from the
last `p` drawn factors with an innovation at every step, and the observable
variables read off the loadings, both error terms drawn. Only the horizon has to
be supplied.

* Added `add_posterior_loglik.dfmodel()`, the pointwise log-likelihood laid out
for WAIC and PSIS-LOO. It is the measurement density at the drawn factors, so it
is the log-likelihood *conditional* on the factor path rather than the marginal
one; the distinction matters for what an information criterion computed from it
means.

* `draw_posterior.dfmodel()` is replaced by `add_posterior_coefficients.dfmodel()`,
following bvartools, which has retired the `draw_posterior` generic.

* The `dfm` class and its methods -- `dfm()`, `plot.dfm`, `summary.dfm`,
`thin.dfm`, `predict.dfm` -- are removed. Posterior draws stay on the `dfmodel`
object, in `object$posterior`, as `coda::mcmc` matrices with one row per draw.

* `Matrix` is no longer a dependency; it was needed only by the R implementation
of the sampler. `bvartools (>= 1.0.0)` is, for the generics
`add_priors`, `add_initial_values`, `add_posterior_coefficients`,
`add_posterior_forecasts` and `add_posterior_loglik`.

* The license is now `GPL (>= 2)` rather than `GPL (>= 3)`, matching bvartools.
Part of what this package now contains was released there under GPL (>= 2), and
that grant cannot be withdrawn, so labelling it GPL (>= 3) here would have given
the same code two different notices without restricting anything. Every
dependency -- bvartools, Rcpp, RcppArmadillo, coda -- is GPL (>= 2) as well.
`inst/COPYRIGHTS` records the BSD-3-Clause portion, which is compatible with
either version, and notes that `bem_dfmdata` is third-party data whose terms are
not this package's to grant.
