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
