# dfmtools (development version)

* **Generalised impulse responses and forecast error variance
  decompositions** for class `favarmodel`. `irf()` gained a `type` argument --
  `"oir"` (the Cholesky default), `"gir"` for the generalised response of
  Pesaran and Shin, which needs no ordering, and `"feir"` for the reduced-form
  response -- and `fevd()` is new, on the `bvartools` generic, so its result
  plots with that package's method.

  A decomposition of a panel series carries one column per state shock plus an
  `"idiosyncratic"` column: the series' own error is white noise, so it enters
  the forecast error variance once at every horizon rather than accumulating,
  but it enters it all the same, and the share of the forecast error the common
  component does not explain is usually the more informative number. An observed
  variable has no such column, being part of the state and measured without
  error. Under `"oir"` the shares sum to one; under `"gir"` they overlap, and
  `normalise_gir = TRUE` rescales the state shares to fill exactly the share the
  state accounts for, leaving the idiosyncratic column alone.

  Note that `"gir"` does not respect the slow-moving restriction: a generalised
  shock to the policy rate moves a slow series on impact, through its loadings
  on the factors, which is what `add_priors(slow = )` was imposed to rule out.
  The two identify different shocks, and `irf()` documents the difference rather
  than picking for you.

* **A vignette, "Identifying a monetary policy shock"**, working the Bernanke,
  Boivin and Eliasz design end to end on `bem_dfmdata`: ordering the panel so
  that slow-moving series identify the factors, `add_priors(slow = )` for the
  contemporaneous restriction, and `irf()` with the funds rate ordered last.

  It is their design rather than their numbers -- theirs is a monthly panel of
  120 series through 2001:8, this is quarterly FRED-QD -- and the vignette says
  where the two part company. The medians reproduce their qualitative pattern,
  including the absence of a price puzzle; the credible bands cover zero at every
  horizon, which the vignette reports rather than buries.

  `knitr` and `rmarkdown` join Suggests, and DESCRIPTION gains a
  `VignetteBuilder` field.

* **`add_priors()` on class `favarmodel` gained a `slow` argument**, the panel
  series assumed not to react to the observed block within the period. Their
  loadings on the observed columns of `lambda` are pinned at zero, which is the
  restriction Bernanke, Boivin and Eliasz identify a monetary policy shock with:
  once the factors are slow-moving, a recursive ordering with the policy rate
  last says something the data has not already been asked to say.

  Series are given by name or by index, and the identifying first `n` may be
  listed without effect -- the identification has zeroed their observed columns
  already, which is also why those `n` should themselves be slow-moving. The
  restriction is a prior rather than a hard zero; `slow = list(series = ...,
  vinv = ...)` sets the precision it is held at, `1e12` by default, which leaves
  a loading at zero to eight decimal places.

  This was previously possible only by writing into `priors$lambda$vinv` at
  hand-computed offsets, which the argument now does.

* **`irf()`** on class `favarmodel`, the response of an observable to a
  recursively identified shock to one element of the state. A shock is an
  element of `s_t = (f_t', y_t')'` -- a factor or an observed variable -- because
  that is where the model's dynamics are; the panel's own error is idiosyncratic
  and propagates nothing. A response is any observable: the `k` panel series
  followed by the `n_obs` observed ones, the order `add_posterior_forecasts()`
  already returns. A panel response is carried through the loadings, which is
  what lets a factor's shock be read on a series that never enters the
  transition.

  Identification is the Cholesky factor of `Q` in the order given by argument
  `order`, which defaults to the state's own, factors before observed variables.
  Nothing in the estimated model tests that ordering -- `Q` is unrestricted by
  construction, which is what the model is estimated for -- so a different one
  gives a different answer from the same draws.

  The return value carries class `bvarirf`, so `plot()` from **bvartools**
  works on it, and `keep_draws = TRUE` returns the posterior itself rather than
  its quantiles.

* **`create_favarmodel()`**, a factor augmented VAR. **Draws are unchanged** for
  the four dynamic factor models: this adds a fifth sampler beside them and
  touches none of their code paths.

  ```
  x_t = lambda_f f_t + lambda_y y_t + e_t,   e_t ~ N(0, R),  R diagonal,
  s_t = sum_j Phi_j s_{t-j} + v_t,           v_t ~ N(0, Q),  s_t = (f_t', y_t')',
  ```

  after Bernanke, Boivin and Eliasz (2005). The observed block `y_t` is part of
  the *state*, not a set of regressors: it appears on the left of the transition
  as well as the right, and `Q` is unrestricted so that its cross block -- the
  correlation between the factor innovations and the shock to the observed
  variables -- can be read. That is what the model is estimated for.

  Five methods: `add_priors()`, `add_initial_values()`,
  `add_posterior_coefficients()`, `add_posterior_forecasts()` and
  `add_posterior_loglik()`, all on class `favarmodel`.

  This is what the previous entry described as arriving without a binding. The
  sampler is now vendored -- `tools/update-bayests-core.R` gained two entry
  points and the closure pulled in `favar_support.h` and `wishart.{h,cpp}` on its
  own, twenty-seven files to thirty-two -- and `src/FavarNormalWishart.cpp` is
  the binding it was waiting for.

  **Three things differ from the dynamic factor models and are worth knowing.**

  *The identification is an identity block, not a unit lower triangle.* The
  leading `n x n` block of the factor loadings is the identity and the observed
  columns of those rows are zero, so the first `n` panel series are the factors
  plus idiosyncratic noise exactly. A dynamic factor model can use a unit lower
  triangle because its `V` is diagonal and the two restrictions together admit
  only the identity rotation; a FAVAR has no diagonal `Q` to offer, so the
  triangle alone would leave the loadings free to wander along a ridge. Order the
  panel so that its first `n` columns are the series you will define the factors
  by.

  *The prior on `Q` is a Wishart*, not a set of independent gammas -- the only
  error block in this package that is not a diagonal. `df` defaults to the width
  of the state and `scale` to the identity of that size.

  *The forecast is wider than the panel.* Each horizon carries the `k` panel
  series followed by the `n_obs` observed ones, `h * (k + n_obs)` columns in all,
  because past the end of the sample the observed variables are no longer data
  and are forecast alongside the factors.

  One thing to watch when reading loadings: with the default `normalize_x = TRUE`
  each panel column is divided by its own standard deviation, so a loading comes
  back on that standardised scale. It is correct there and not comparable with
  one implied by the raw data until the normalisation is undone.

  Bernanke, B. S., Boivin, J., & Eliasz, P. (2005). Measuring the effects of
  monetary policy: A factor-augmented vector autoregressive (FAVAR) approach.
  *The Quarterly Journal of Economics, 120*(1), 387-422.

* **Vendored BayesTS core refreshed.** **Draws are unchanged**, for all four
  dynamic factor models. Upstream's own fingerprint comparison reports 76
  fixtures unchanged and none moved over the shared code this picks up.

    Two things arrive that nothing here calls yet, both from upstream's work on a
    factor augmented VAR. `core/algorithms/chan_jeliazkov_2009.cpp` gains a
    second entry point, `chan_jeliazkov_2009_conditional`, which holds the
    trailing elements of every state column at observed values instead of drawing
    them -- what a state vector that is part data needs -- and
    `core/models/dfm_support.h` gains `draw_conditional_factor_path()` beside
    `draw_factor_path()`, on the same one-block-or-stack contract and with the
    same shift convention. The four samplers here reach neither. Their arrival
    split the band sampler into an assembly, a conditioning step and a draw;
    upstream verified that split moved no fingerprint.

    The rest is type surface for a model this package does not yet have.
    `FavarNormalWishart` is a factor model and so belongs here rather than in
    bvartools, which skips it; it is *not* vendored, because there is no binding
    and no R entry point, and compiling a sampler nothing can call is the trade
    `src/core/VENDORED.md` declines. Adding it later is two lines in that
    script's `entry` list. What is compiled in meanwhile is its `Input`,
    `Initial` and `Draws` structs, its `validate()`, `n_obs_factors` with
    `n_state()`, `n_favar_lambda()` and `n_favar_a()`, and `TrainData::f_obs` --
    all in files every model shares and copied whole.

* **Both drifts at once**, with `create_dfmodel(error = "sv", tvp = TRUE)`. The
two arguments are now independent, so the four combinations select the four
algorithms `DfmNormalGamma`, `DfmNormalStochvol`, `DfmTvpGamma` and the new
`DfmTvpStochvol`. In that last one the measurement equation is
$x_t = \lambda_t f_t + u_t$ with $u_t \sim N(0, U_t)$ and the transition
$f_t = \sum_i A_{i,t} f_{t-i} + v_t$ with $v_t \sim N(0, V_t)$: nothing in it is
held fixed but the normalisation.

    Carrying both is not redundant, which is the point of having it. A series
    whose loading fell looks like a series whose idiosyncratic variance rose, and
    a period of common turbulence looks like a transition that changed; a model
    with only one of the two drifts has to explain the other with what it has.
    The two are separated because the loading paths are drawn weighted by their
    own series' volatility period by period — the periods in which a series was
    quiet identify its loading path and the periods in which it was wild largely
    do not, while the path is free to move between them.

    Nothing new to learn at the interface: the specification is the union of the
    two that already existed. `add_priors()` takes the state equation
    (`shape`, `rate` on top of `vinv`) for `lambda` and `a` and the six-element
    stochastic volatility specification for `u` and `v`, in one call. All four
    groups then carry a `shape` and a `rate`, at four different widths, and each
    pair belongs to the random walk of the block it sits in.
    `add_initial_values()` returns four paths. In the posterior, `lambda` and `a`
    are $t$ times wider because the coefficients drift and `u_sigma_inv` and
    `v_sigma_inv` are $t$ times wider because the variances do — this is the only
    model here in which all four are paths. `add_posterior_forecasts()` holds
    every one of them at its last in-sample period;
    `add_posterior_loglik()` scores every period under its own $\lambda_t$ *and*
    its own $U_t$.

    **Draws are unchanged** for the three models that existed before. The
    vendored core was refreshed for this and the change it brought to shared code
    is one merge — the factor path draw, which was three near-copies of a period
    shift convention and is now one function serving all four models. Verified
    upstream with the fingerprint comparison rather than argued: 74 fixtures
    unchanged, none moved, over the merge and again over the new sampler.

    One thing that changed for a caller: `create_dfmodel(error = "sv",
    tvp = TRUE)` used to be an error naming the combination as unavailable, and
    is now a model.


* **Time varying loadings and transition coefficients**, with
`create_dfmodel(tvp = TRUE)`. The measurement equation becomes
$x_t = \lambda_t f_t + u_t$ and the transition
$f_t = \sum_i A_{i,t} f_{t-i} + v_t$, with every freely estimated element of
$\lambda_t$ and every element of the $A_{i,t}$ a random walk of its own, drawn as
a whole path with the simulation smoother of Durbin and Koopman (2002). The new
algorithm is `DfmTvpGamma`, beside `DfmNormalGamma` and `DfmNormalStochvol`, and
`tvp` names the coefficients exactly as it does in `bvartools`' VAR and VEC
models -- it moves *both* blocks, and a model in which only the loadings drifted
would be a different one.

    What this is for is the assumption a factor model makes most often and
    defends least: that a series' exposure to the common component held over the
    whole sample. A series can enter or leave that component -- a sector
    reorganised, a country's trade opening -- without anything about the factor
    itself changing, and a constant-loading model has nowhere to put it except
    the idiosyncratic variance, which then carries it as noise the series is
    credited with throughout, including in the periods where the exposure did
    hold. Drift in the transition is the other half: the persistence of the
    common component is what a forecast from it runs on, and it is not a constant
    of nature either.

    The leading $N \times N$ block of $\lambda_t$ does not move, because it is
    not estimated. Only the product $\lambda_t f_t$ is identified, so that block
    is the normalisation that pins the rotation and the scale of the factors;
    letting it drift would let both wander over the sample, and a loading path
    would then describe the normalisation as much as the exposure it is read as.

    What changes for a caller of the existing models: nothing. `tvp = FALSE` is
    the default and the draws of `DfmNormalGamma` and `DfmNormalStochvol` are
    unchanged -- the only shared code the addition touched is one shape dispatch
    in the vendored core, verified against the upstream fingerprint comparison as
    72 fixtures unchanged and none moved. For `tvp = TRUE`, arguments `lambda`
    and `a` of `add_priors()` take `shape` and `rate` on top of `vinv`, which is
    the inverse gamma on the variance of the state innovations; `vinv` itself
    then describes the state of the period before the sample rather than a
    coefficient that holds throughout. There are no defaults for the pair, so a
    call that forgets is told which elements are missing. `add_initial_values()`
    returns `lambda` and `a` as matrices with one column per period, beside
    `lambda_init`/`a_init` and `lambda_sigma_inv`/`a_sigma_inv`. In the posterior,
    `lambda` and `a` widen from one number per coefficient to a whole path --
    $M N \times T$ and $N^2 p \times T$ per draw, the periods stacked within a
    row -- and each gains a `sigma` element holding the variance of its state
    innovations. `add_posterior_forecasts()` holds both at their last in-sample
    period, as every time varying model in this family does, while
    `add_posterior_loglik()` scores every period under its own $\lambda_t$.

    Only `error = "gamma"` for now. `error = "sv"` with `tvp = TRUE` is refused
    rather than silently estimated as something else.


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
