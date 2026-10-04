# dfmtools (development version)

* **Deterministic terms in every factor model.**
  - **The arguments.** `create_dfmodel()` and `create_favarmodel()` take
    `deterministic` (`"none"`, the default, `"const"`, `"trend"` or `"both"`)
    and `seasonal`, as bvartools' `create_bvarmodel()` does, and build the same
    terms: `const`, `trend` and `season.1` onwards, in `data$deterministic`.
    Until now a panel had to be demeaned and deseasonalised before it was
    passed.
  - **Where they enter.** In the measurement equation,
    `x_t = lambda f_t + C d_t + u_t`, so the factors are deviations from them.
    A FAVAR's observed variables deviate from `C_obs d_t` as well, and its
    transition runs over the deviations. `C` does not drift under `tvp = TRUE`.
  - **Priors, starting values and draws.** `add_priors()` takes `c` (and
    `c_obs` for a FAVAR), `add_initial_values()` starts both at their least
    squares values, and the draws are `posterior$c$coeffs` and
    `posterior$c_obs$coeffs`, named `"<series>.<term>"`. The constant and the
    trend trade off against the level of a persistent factor and mix slowly, so
    such a model needs a longer burn-in.
  - **The core.** The vendored BayesTS core moves to `065ed77`, which adds the
    terms to all five samplers. The draws of a model without deterministic
    terms are unchanged.

* **Forecasts take bvartools' route.**
  - **The steps.** `add_forecast_input()` sets the horizon and the
    deterministic terms of the forecast periods, continuing them from the
    sample or taking them as `deterministic = `; `prepare_forecast_input()`
    returns them without adding them. `add_posterior_forecasts()` then
    simulates, and the new `predict()` methods collect the draws as a
    `'bvarprd'`, in the data's units, which bvartools' `plot()` draws.
  - **Forecast errors.** `add_forecast_errors()` and `get_forecast_errors()`
    have methods for both classes.
  - **Deprecated.** `add_posterior_forecasts(model, n_ahead = h)`, and calling
    it without a horizon, still work through `add_forecast_input()`, and warn.
    The forecasts are the same.

# dfmtools 0.1.0

* **Factor models can be written to BayesTS model files and read back.**
  `write_to_hdf5()` has methods for a `'dfmodel'` and a `'favarmodel'`. The
  file can be estimated, forecast and scored by the BayesTS command line, with
  bvartools' `bayests_files()` for instance. bvartools'
  `read_model_from_hdf5()` returns the model it was written from. A
  `'modellist'` holding factor models beside VARs and VECs is written and read
  back with bvartools' folder functions, in its order and with every class
  intact.
  - **Translation.** The file follows BayesTS's layout, not the R object's:
    - `m` and `n` become `k` and `n_factors`. BayesTS would read `m` and `n`
      as a VAR's exogenous variables and deterministic terms.
    - The panel becomes `/data/train/y`. The deterministic terms become
      `/data/train/x`, with their number as `/model/n`, and the terms over
      the horizon become `/data/forecast/x`.
    - `u` and `v` become `/priors/u_sigma` and `/priors/v_sigma`.
    - The free loadings are put in the sampler's order, row by row. This uses
      the samplers' own permutation, so a model read back gives the same draws
      as the original, under a prior that is not the same in both orders.
  - **Normalisation.** The panel's centres and scales are kept, so a model
    read back can still be interpreted on the scale of its data.
  - **Resampled draws.** The draws that `add_sign_zero_restrictions()` keeps
    no longer gain coda's placeholder column names `var1`, `var2`, ….
  - **Dependency.** This needs a bvartools that exports
    `write_bayests_tree()` and `from_bayests_tree()`.

* **Factor augmented VARs can be identified by sign and zero restrictions.**
  - **The function.** `add_sign_zero_restrictions()` has a method for a
    `'favarmodel'`. It rotates the shocks of the state equation with the
    algorithm of Arias, Rubio-Ramírez and Waggoner (2018), through
    bvartools' `arias_rubio_ramirez_waggoner_2018()`, and resamples the
    posterior by the importance weights. It takes `max_tries`, as
    bvartools' method does.
  - **Using the rotations.** `irf()` and `fevd()` read the rotations under
    `type = "sign"`.
  - **What can be restricted.** Restrictions are on elements of the state:
    the factors, each named after the panel series that identifies it, and
    the observed variables. The other panel series respond through
    draw-specific loadings, which the algorithm cannot restrict, and are
    refused.
  - **Dependency.** dfmtools now needs bvartools 1.0.0.

* The vendored BayesTS core moves to `696dc7d`, which adds mixed-frequency and
  missing-data estimation, also beside i.i.d. variables, new prior options,
  among them the normal-gamma prior with a fixed or a drawn theta, and the
  structural quantile VAR to the VAR models. None of it
  applies to a factor model, and the draws of every model are unchanged.

* **`expected_model_size()` says how large a factor model will be before it is
  estimated.** The generic belongs to bvartools, and the methods for a
  'dfmodel' and a 'favarmodel' size every block of the posterior from the
  specification alone, as the methods there do for a VAR. The factors are
  drawn for every period, and time-varying parameters and stochastic
  volatility multiply the loadings, the transition and the error precisions by
  the number of periods as well: a one-factor model of `bem_dfmdata` with
  drifting loadings needs 1.8 GB at 5000 draws. `add_posterior_coefficients()`
  warns before it starts when the model will exceed
  `options(bvartools.size_warning)`, as it does for the models of bvartools.

* `add_priors()` takes the non-centred prior `omega_v` for every random walk a
  model with `tvp = TRUE` or `error = "sv"` has: the loadings, the factor
  transition and the log-volatilities of both error terms. It replaces `shape`
  and `rate` for the block it is given on, and puts a normal prior on the signed
  standard deviation of the state innovations rather than an inverse gamma on
  their variance, following Frühwirth-Schnatter and Wagner (2010). The four
  blocks are independent, so a model may take it for some and keep the gamma
  prior for the others.

* New method `time_variation_test()` for `dfmodel` objects, which reports the
  Savage-Dickey Bayes factors of Chan (2018) for time variation, one per state
  and one per block, for whichever blocks were estimated under `omega_v`. It is
  the `bvartools` generic, so the tables of a factor model and of a VAR read
  alike.

* The vendored BayesTS core is refreshed to upstream `d5ca82c`, which is where
  the two above come from. The sampler files of `DfmTvpGamma` and
  `DfmTvpStochvol` moved for the new parameterisation and a new file,
  `noncentred_support.h`, arrived with it; the core also gained a refusal for
  `/model/n_iid`, a VAR restriction a factor model does not offer. *Draws are
  unchanged* for a model that does not set `omega_v`: upstream fingerprints all
  120 of its existing fixtures identically across the commit.

* The vendored BayesTS core is refreshed to upstream `d92e581`. Its samplers
  now refuse a NaN or an infinity in the data, the priors, the initial values
  and the realised values a forecast is scored against, and name the input
  that holds it. Before, such a value failed later, deep in the numerics, or
  in a forecast or a score was not caught and came back as NaN. An NA in the
  series still stops `add_initial_values()` before the sampler is reached,
  as it did before. The other files that moved change only comments.
  *Draws are unchanged*: upstream fingerprints every fixture identically
  across the commit.

* The vendored BayesTS core is refreshed to upstream `f8b42a1`. Only
  `inputs.cpp` moved, and only a comment in it: the description of
  `validate_tvp_block()` is back above that function, from where an addition
  had separated it. *Draws are unchanged*.

* The vendored BayesTS core is refreshed to upstream `94f81de`. Three files
  moved. `kalman_durbin_koopman_2002.cpp` skips the products with the
  identity transition of a random walk, which makes the smoother faster and
  is exact. `forecast_states.h` guards the square root forecast errors are
  drawn through, which could give NaN from a badly conditioned precision in
  the VAR and VEC models, and `favar_normal_wishart.cpp` now calls it in place
  of its own copy. That copy already guarded against the same thing, and the
  two differ only in what a negative eigenvalue of Q becomes -- zero now, its
  absolute value before -- which is a rounding error either way. *Draws are
  unchanged*: none of upstream's factor model fixtures fingerprints
  differently across either commit.

* The vendored BayesTS core is refreshed to upstream `cea124b`, which adds the
  non-centred prior for `VecTvpGamma`: two members of that model's draw struct
  in `results.h` and a comment in `priors.h`, neither read by a factor model.
  *Draws are unchanged*; upstream fingerprints every fixture identically across
  the commit.

* The vendored BayesTS core is refreshed to upstream `ffbd208`. Four files
  moved -- `inputs.h`, `priors.h`, `results.h` and `inputs.cpp` -- for the
  non-centred random walks of `VarTvpGamma` and `VecTvpStochvol` and for the
  VEC models' cointegration prior, none of which a factor model uses. *Draws
  are unchanged*: every factor model fixture fingerprints identically across
  each of those upstream commits.

* The vendored BayesTS core is refreshed to upstream `4a64082`. Nothing a
  factor model draws changes: the seven files that moved carry a warning and
  three refusals for the variable selection the VAR and VEC samplers offer,
  which no factor model does, and the building blocks of a non-centred
  stochastic volatility that only `VarTvpStochvol` uses -- additions beside
  the mixture draw the factor models call, not changes to it. *Draws are
  unchanged*: upstream's fingerprint recording, which includes every factor
  model fixture, is identical before and after each of those commits.

* **A univariate block passed as a plain vector broke a factor augmented VAR.**
  `create_favarmodel()` left `x` and `y` as they arrived, so a `ts` vector --
  which is what a single observed variable is, unless a caller thought to keep
  it a one-column matrix -- reached the binding as a vector and came back as
  Rcpp's `Not a matrix.`, which `add_posterior_coefficients()` turned into
  `error = TRUE`. All the caller then saw was that there were no draws.
  `create_dfmodel()` has always converted; both now use the same helper, which
  keeps the `tsp` and names columns that R never named.

* `Depends` requires `R (>= 4.0.0)` rather than `(>= 3.5)`. `src/Makevars` sets
  `CXX_STD = CXX17` and 3.5 predates a toolchain that reliably honours it, so
  the old floor was a claim nothing tested -- the check matrix goes back only
  to `oldrel-1`.

* `Date` in `DESCRIPTION` and `date-released` in `CITATION.cff` are current and
  agree, as `.github/CONTRIBUTING.md` asks of them, and `.Rbuildignore` no
  longer names a `CONTRIBUTING.md` that has not been at the top level for some
  time.

* **A test sample with a missing period was scored against the wrong periods.**
  `add_predictive_loglik()` dropped incomplete rows with `na.omit()`, which on a
  time series is `na.omit.ts()`: where it drops leading `NA`s it advances the
  `tsp` of what it returns, and the period labels were computed from the
  original start, so every realised value was attributed one period early.
  Interior `NA`s it refused outright, and a wholly missing sample was an error
  rather than no periods to score. Each row is now matched to the period it came
  from before anything is dropped, and what is scored is the unbroken run from
  the first forecast period -- a gap ends the run rather than being skipped,
  because the filter conditions each period on the one before it.

* `add_posterior_coefficients()` takes `verbose`. The sampler has been able to
  report its progress since the C++ core arrived -- a percentage, at most one
  line per percent -- and no argument reached it, so a chain of the default
  length said nothing for minutes. It is still silent by default, and reporting
  changes no draw.

* The documentation of the seed no longer lets the reproducibility of
  `add_posterior_coefficients()` be read as covering the steps after it.
  `add_posterior_forecasts()` is not seeded: its draws come from R's generator
  as it stands, so on a cluster they depend on which worker took the model, and
  the same list forecast on one worker and on four does not give the same paths.

* Tests for what was never run: all four samplers are scored rather than the two
  with constant error variances, which is the code the one defect this filter
  has had was in; the factor augmented VAR's methods have their refusals
  covered; the shipped `inst/COPYRIGHTS` is checked against the files actually
  vendored, which only `tools/update-bayests-core.R` did and only with an
  upstream checkout to hand; and the package's unload hook runs in a process of
  its own. Coverage of `R/` goes from 91.9% to 95.1%.

* `data-raw/bem_dfmdata.R` says what it does and checks the one thing it cannot
  see: that the headerless panel and the FRED-QD file it takes the series names
  from describe the same columns in the same order. A mismatch would have named
  all 196 series wrongly in silence. The script reproduces the shipped object
  exactly.

* **The prior arguments are checked for what they are, not only for where they
  lie.** `spec$shape < 0` on the character `"a"` compares as strings and is
  `FALSE`, so every range check in this package passed exactly the argument it
  was written to catch. `add_priors()` on class `favarmodel` took a character
  `lambda$vinv` without a word -- `diag()` coerced it to `NA` and only warned,
  leaving an all-`NA` precision matrix -- and read `v$df` with the same string
  comparison; the stochastic volatility fields of class `dfmodel` took one too,
  and the complaint arrived from `add_initial_values()` as `missing value where
  TRUE/FALSE needed`. Every prior field of both classes now goes through one
  helper that checks the type first and names the argument, and the two
  `add_priors()` methods share it, along with the builder of the idiosyncratic
  gamma priors, so they cannot drift apart again.

* **A flat prior is refused where it is written.** Every precision, and every
  gamma a starting value is drawn from, must now be larger than 0 rather than
  at least 0 as the documentation used to say. `vinv = 0` was accepted and then
  had no variance to draw a starting value with: a dynamic factor model stopped
  inside `chol()` with "the leading minor of order 1 is not positive", and a
  factor augmented VAR drew `rnorm(sd = Inf)`, put `NA` starting values in front
  of the sampler and came back complaining about `NaN` from
  `chan_jeliazkov_2009`. A gamma of shape 0 was the same story one function
  later. The one place a zero is still a number rather than a mistake is the
  loading block of a model with a single series and a single factor, which has
  no freely estimated element for the prior to be about, and `shape` of a
  stochastic volatility specification, which nothing draws from.

* `v$scale` of a factor augmented VAR must be positive definite, not only
  symmetric. A matrix of zeros passed the symmetry test trivially --
  `max(abs(0 - 0)) <= 0` -- so the one scale that is certainly wrong was the
  one that got through.

* A wrong-sized loading prior or starting value is reported by the size it
  actually has. The bindings set these only when the size already matched, so
  the sampler's validator described what it had been given rather than what the
  caller wrote: `prior precision of lambda must be 7x7, got 0x0` for a 3 x 3.

* `add_posterior_coefficients()` keeps what a caller has put on the model
  object. It was rebuilt from `data`, `model`, `initial`, `priors` and
  `posterior`, which silently dropped everything else; it is now assigned into,
  and the R method clears a stale `posterior` and `error` before the call.

* `n_ahead` of `add_posterior_forecasts()`, `irf()` and `fevd()` must be a whole
  number. `n_ahead = 2.7` was accepted and quietly behaved as 2.

* `release-version.yaml` no longer truncates a tag at the first dash. An R
  version number may use `-` as a separator, so `v0.2-5` would have been checked
  as `0.2`; only a recognised pre-release suffix is stripped now.

* **Starting values of the error precisions came from the wrong distribution.**
  `add_initial_values()` on class `dfmodel` drew each element of `uinv` and
  `vinv` as the reciprocal of a gamma drawn with the rate inverted, rather than
  from the Gamma(`shape`, `rate`) prior that `add_priors()` had stored on the
  precision and the sampler draws it from. At the documented defaults, shape 5
  and rate 4, that is a mean of 1/16 where the prior mean is 5/4 -- a factor of
  twenty, and an initial idiosyncratic *variance* near 16 on a panel the
  normalisation has already set to variance one. Only the starting values were
  affected, so this cost burn-in rather than correctness, and a chain long
  enough to have converged gives the same posterior as before. Short chains move.
  `add_initial_values()` on class `favarmodel` and the error precisions of a
  bvartools VAR were already drawn the right way; all three now go through one
  helper, which also refuses a `shape` of zero rather than starting from the
  singular precision `rgamma()` returns for it.

* **`irf()` on class `favarmodel` returned a transposed object for
  `n_ahead = 0, cumulative = TRUE`.** Accumulating over a single horizon left
  `apply()` with nothing to keep the horizon dimension for, so a draws-by-one
  matrix came back as one-by-draws and the credible interval was then taken
  over the draws rather than over the horizon. `fevd()` already guarded the
  same case where it accumulates.

* **A test sample is put on the model's scale before it is scored.**
  `create_dfmodel()` normalises the panel by default, so a forecast and the
  density `add_predictive_loglik()` computes are on the standardised scale --
  but `test_sample` was passed through in the data's own units and scored
  against it as it stood, silently and with nothing in the result to say so.
  Realised values are now centred and scaled with the *estimation sample's*
  moments, and `data$test$x` holds what was scored, on that scale. A model
  created with `normalize_x = FALSE` is unaffected and is scored in the data's
  units, as it was estimated in them; the two differ by the Jacobian of the
  normalisation, so densities are comparable across models only when the models
  normalise alike. A new section of `?add_predictive_loglik.dfmodel` says this.
  Scores computed with an earlier version on a normalised model are wrong and
  worth recomputing.

* `posterior$loglik` is a `coda::mcmc` object rather than a bare matrix, for
  both classes, so a thinned model's pointwise log-likelihood carries the same
  draw labels as its coefficients -- which every other block of draws in the
  package already did. `waic()` and `loo()` take it as the matrix it still is.

* Of the steps after estimation, `add_posterior_forecasts()`,
  `add_posterior_loglik()` and `add_predictive_loglik()` now report a failure
  and stop, for both classes, rather than three of them returning the model
  with `error = TRUE`. They are cheap and derived, and what goes wrong in them
  is nearly always the call rather than the chain, so a message beats a flag.
  `add_posterior_coefficients()` is unchanged and still carries on, which is
  what makes a list of models lose only the one that failed.

* `create_dfmodel()` and `create_favarmodel()` check `p`, `n`, `iterations` and
  `burnin` where the caller names them, as they have always checked `thin`. A
  non-integer or non-numeric `p` used to travel as far as a matrix
  multiplication -- `"a" < 0` compares as strings and is `FALSE` -- and a
  non-positive `iterations` as far as the sampler, a priors and starting values
  round trip later.

* `add_initial_values()` on class `dfmodel` says when `priors` is missing
  instead of failing inside `nrow(NULL)`, as the `favarmodel` method already
  did. `add_priors()` on the same class says which argument a `NULL` block was,
  instead of failing inside `diag()`, and reports a short `u` or `v` by the
  field it is short of rather than by its length.

* Documentation fixes: `README.md` no longer contains the four R errors a
  render against a stale installed package baked into it; the
  `tvp-and-stochastic-volatility` vignette reads the forecast draws from
  `posterior$forecast$forecasts` rather than from the `posterior$forecast` that
  stopped being a matrix; `create_dfmodel()` names the returned data element
  `x`, which is what it is; two help pages show bold text rather than the
  literal asterisks around it; the three `favarmodel` methods that had no
  examples have them; and the places that still asked for `bvartools (>= 1.0.0)`
  now agree with the `DESCRIPTION` that has required `(>= 0.3.0.9000)` since
  September. `tools/render-readme.R` and `vignettes/precompile.R` both stop on a
  failed chunk and render in English, and `.github/CONTRIBUTING.md` says to
  regenerate both from the installed package.

* A `test-coverage.yaml` workflow, the same one bvartools runs.

* **Vendored BayesTS core refreshed: a scored forecast keeps its filter positive
  semi-definite.** `add_predictive_loglik()` on class `dfmodel` could stop with
  "the one step ahead forecast variance of a scored period is not positive
  definite" on input there was nothing wrong with, on some toolchains and not
  others. The filter the score runs updated its state covariance in the short
  form `P - K F K'`, a difference of two positive semi-definite matrices. The
  innovation variance of a transition of order above one is singular -- it
  carries the factor variances in its leading block and zeros elsewhere -- so
  the lagged blocks of that covariance are never refreshed and rounding error
  accumulates in them instead of being flooded out, until the matrix drifts
  indefinite and the check rejects a forecast variance that is mathematically
  fine. Whether the drift crosses zero depends on which BLAS rounds which way,
  which is why upstream's Windows job found it and its two others did not. The
  update is now the Joseph form, `(I - K Z) P (I - K Z)' + K R K'`, a sum of two
  positive semi-definite terms, and the density is evaluated through a Cholesky
  of the forecast variance rather than an LU determinant and a general solve.

    **Scores change by a rounding error** and nothing else does. The four
    dynamic factor samplers were fitted from a pinned seed against the package
    built before and after the refresh, through the coefficients, the pointwise
    log likelihood, a three-step forecast and the score of it: every
    `posterior$forecast$loglik` moved in its last bit, by at most 2.9e-16
    relative, and every coefficient draw, forecast path and pointwise log
    likelihood is bit-identical. The factor augmented VAR is untouched, having
    no score to compute.

    The refresh takes the core to upstream `7d7c6ca`, which is BayesTS 0.3.0
    and that fix. What 0.3.0 itself adds is a pair of discounted time varying
    VAR and VEC models, which are not vendored here and reach this package only
    as declarations in `bayests/inputs.h`, `priors.h`, `results.h` and `spec.h`
    that no sampler here reads. The vendored set is still the same 35 files, so
    `inst/COPYRIGHTS` lists the same names -- though `src/core/VENDORED.md` had
    gone on calling the copy BayesTS 0.2.0 through the two refreshes before this
    one, and now records the commit it is actually at and the three files those
    refreshes added.

* **A forecast can be scored against what its horizon realised.**
  `add_predictive_loglik()` takes a fitted model with a forecast and a test
  sample, and adds `posterior$forecast$loglik`: one row per draw and one column
  per scored period, each column the log density of that period's realised
  observation given the ones before it. The log of the mean of a column's draws
  is the one step ahead predictive density, and those sum over the periods to
  the log predictive likelihood of the realised stretch. What the model was
  scored against is kept in `data$test$x`, so the same call can be made again
  without the sample.

  A factor model is scored by **filtering**, which is not how bvartools scores a
  VAR. A VAR reaches its realised history through its regressors, so its score
  is its own pointwise log-likelihood on another sample. A factor model's
  history reaches the density through the latent factors, and its in-sample
  log-likelihood conditions on the factors the sampler drew -- of which there
  are none outside the sample. So at every scored period the realised
  observation updates the distribution of the factors before the next is
  predicted, and the column is that period's prediction error decomposition.
  Where the loadings, the transition or the volatilities drift they take a step
  per scored period, in the order `add_posterior_forecasts()` steps them; only
  the free elements of the loading matrix walk, the identifying block being what
  fixes the rotation and the scale.

* **`posterior$forecast` is a group rather than a matrix of draws.** The
  simulated paths are now at `posterior$forecast$forecasts`, with the score
  beside them at `posterior$forecast$loglik`. **This breaks code that reads
  `posterior$forecast` as a matrix.** It is the layout bvartools uses and the
  one the model file writes as `/posterior/forecast`, and the members are named
  after what they hold because all of them are draws.

* **Lists of models can be simulated on several cores.** A `modellist` that
  `create_dfmodel()` or `create_favarmodel()` returns for more than one lag
  order or number of factors is simulated by the list methods of bvartools,
  and their argument `cores` works for it: `add_posterior_coefficients()`,
  `add_posterior_forecasts()` and `add_posterior_loglik()` share the models out
  over that many worker processes, which load dfmtools for the methods they
  need. Coefficient draws are seeded per model and equal those of a process
  running with one thread, as the workers do; for a factor augmented VAR that
  can differ in the last digits from a session running several. Forecasts
  draw from worker streams that `set.seed()` reproduces. Needs the bvartools
  whose workers load the packages of the models' methods; with earlier
  development versions a list of these models stopped with "no applicable
  method". `tests/testthat/test-parallel.R` runs both kinds of list on two
  workers.

* **Every model carries the seed of its posterior simulation, as in
  bvartools.** `add_initial_values()` on class `dfmodel` or `favarmodel` stores
  a seed, drawn from R's random number generator, as `object$model$seed` unless
  the model has one, and `add_seed()`, the bvartools generic, gains methods for
  both classes that replace it. The internal samplers of
  `add_posterior_coefficients()` draw with that seed: R's generator is set to
  it, with R's default kinds, for the simulation and put back as it was
  afterwards. **This changes the draws of existing code**, not their
  distribution: a `set.seed()` between `add_initial_values()` and
  `add_posterior_coefficients()` no longer has an effect, and earlier results
  do not reproduce draw for draw. `set.seed()` before `add_initial_values()`
  still makes a script reproducible, and a model now draws the same whichever
  worker of a cluster simulates it. A model without a seed draws from R's
  generator as before, and so do the forecasts of `add_posterior_forecasts()`.

    `add_seed()` on a list returned by `create_dfmodel()` or
    `create_favarmodel()` gives its models the seeds `seed`, `seed + 1`, ...,
    and `add_initial_values()` on such a list already gives each its own. The
    `add_seed()` generic needs a bvartools that has it, and seeding the models
    of a list one whose list method counts the models of other packages.

* **Vendored BayesTS core refreshed: forecasts carry the drift forward.**
  `add_posterior_forecasts()` forecast the three dynamic factor models whose
  states drift -- `error = "sv"`, `tvp = TRUE`, and both -- from each draw's
  loadings, transition and volatilities in the last sample period, held for every
  forecast period. That is the forecast of a model whose drift stops where the
  sample does: its intervals left out the drift, and a held volatility
  understated the expected variance of every period after the first. Each draw's
  random walks now take one step per forecast period, by the innovation variances
  the sampler drew for them; only the free loadings move, and the identifying
  block stays fixed. The new argument `forecast_states = "hold"` gives the old
  forecasts and is stored in `model$forecast_states`. Under stochastic volatility
  `add_posterior_coefficients()` now also keeps `sigma` in `u_sigma_inv` and
  `v_sigma_inv`, the variances of the two groups of log-volatility innovations; a
  model fitted with an earlier version lacks them and forecasts only with
  `forecast_states = "hold"` until it is fitted again. **Forecasts change** for
  those three models. Their coefficient draws and log likelihoods, and everything
  about `DfmNormalGamma` and the FAVAR, are unchanged, which upstream verified by
  fingerprinting every sampler before and after. The vendored set grows by one
  file, `core/models/forecast_states.h`, now listed in `inst/COPYRIGHTS`.

* **Vendored BayesTS core refreshed: time varying coefficients leave their
  starting values.** **Draws change for the two time varying dynamic factor
  models**, `tvp = TRUE` with either error specification, and for nothing else.
  Each of the five samplers was fitted from a pinned seed against the package
  built before and after the refresh, through the coefficients, the log
  likelihood and a six-step forecast. Every posterior block of the two constant
  dynamic factor models and of the FAVAR is bit-identical; every block of the two
  time varying ones moved, and all of it is finite.

    The time varying samplers drew each coefficient path with the simulation
    smoother centred on the previous draw of the state before the sample, and
    with the random walk's own innovation variance as the prior covariance of the
    first period. They drew that state only after the variance. Each step was a
    valid conditional, but together they tied the path's first period to the
    state before it with that variance, and a coefficient meant to drift slowly
    has a small one: the chain could return its starting path as the posterior.
    Upstream now integrates that state out of the first period, drawing the path
    under N(mu_0, V_0 + Sigma), and draws the state given the path before the
    variance conditions on it. Chains started far apart now agree. The change is
    in `core/models/dfm_tvp_gamma.cpp`, `core/models/dfm_tvp_stochvol.cpp` and
    `draw_random_walk_state()` in `core/models/model_support.h`.

    Integrating the state out inverts its prior precision, so `core/inputs.cpp`
    now requires that precision to be positive definite for every time varying
    block. `bayests/priors.h` changes only a comment, on the cointegration space
    prior of the error correction models.

    The vendored set is still the same 32 files, so `inst/COPYRIGHTS` is
    unchanged.

* **`add_priors()` on class `dfmodel` refuses `lambda$vinv` or `a$vinv` of 0 for
  a model with time varying coefficients.** There it is the prior precision of
  the state before the sample, which the refreshed sampler inverts, so a flat
  prior has no variance to give and would be refused when the sampler starts.
  `add_priors()` refuses it first, naming the argument. The check applies only to
  a block with elements, so the loadings of a single series on a single factor,
  which have none, still accept 0. A constant model accepts `vinv = 0` as before.

* **Vendored BayesTS core refreshed: upstream's value checks on the inputs.**
  **Draws are unchanged**, for all four dynamic factor models and the FAVAR,
  verified the same way as the refreshes below. Each of the five was fitted from
  a pinned seed against the package built before and after it, through the
  coefficients, the log likelihood and a six-step forecast. Every posterior block
  of every one of them is bit-identical.

    Upstream's validators used to check the shape of an input and nothing else.
    `core/inputs.cpp` now also refuses values no prior can have, and some of
    those checks are reached from the factor models: a gamma shape or rate that
    is negative or not finite, a log-volatility offset or initial innovation
    variance that is not positive, an initial precision of a time varying block's
    innovations with anything off its diagonal, and a prior precision or Wishart
    scale that is not symmetric. `add_priors()` and `add_initial_values()`
    already refuse or never build every one of these, so no model built from
    their arguments is turned away. A prior precision that is not symmetric
    written into `priors` by hand used to run, and is now refused when the
    sampler starts. `core/models/model_support.h` gains a check on the rows of
    the forecast regressors, which no factor model has, and a comment.

    The vendored set is still the same 32 files, so `inst/COPYRIGHTS` is
    unchanged.

* **`add_priors()` on class `favarmodel` refuses a `v$scale` that is not
  symmetric.** It checked the dimensions only, so a Wishart scale that was not
  symmetric, or had a missing element, ran a chain. The vendored core now refuses
  it when the sampler starts; `add_priors()` refuses it first, naming the
  argument, with the same relative tolerance of `1e-8`. Draws are unchanged for
  every scale it accepts.

* **Vendored BayesTS core refreshed: upstream's BVS, SSVS and TVP log likelihood
  audit.** **Draws are unchanged**, for all four dynamic factor models and the
  FAVAR, verified the same way as the refreshes below. Each of the five was
  fitted from a pinned seed against the package built before and after it,
  through the coefficients, the log likelihood and a six-step forecast. Every
  posterior block of every one of them is bit-identical.

    Nothing this package calls moved. The audit's fixes are in the VAR and VEC
    samplers, none of which are vendored here, and only two vendored files take
    any of it. `core/models/model_support.h` gains `precision_stride()`,
    `half_log_det_precision()` and `stacked_identity()`, and `fill_psi_path()`
    now writes a time-varying Psi as a stack of k x k blocks rather than a block
    diagonal. No factor model calls any of the four. `core/algorithms/stochvol_mixture.h`
    changes only a comment, on how the mixture means are held.

    The vendored set is still the same 32 files, so `inst/COPYRIGHTS` is
    unchanged.

* **Vendored BayesTS core refreshed: cointegration priors must be symmetric.**
  **Draws are unchanged**, for all four dynamic factor models and the FAVAR,
  verified the same way as the refreshes below. Each of the five was fitted from
  a pinned seed against the package built before and after it, through the
  coefficients, the log likelihood and a six-step forecast. Every posterior block
  of every one of them is bit-identical.

    Nothing this package calls moved. Upstream's `core/inputs.cpp` now refuses a
    prior precision of the cointegration space, a prior precision of beta before
    the sample, or a transition `P_tau` that is further from symmetric than
    rounding explains. The check is reached only from the error correction
    validators, none of whose samplers are vendored here.

    The vendored set is still the same 32 files, so `inst/COPYRIGHTS` is
    unchanged.

* **The samplers thin: `create_dfmodel()` and `create_favarmodel()` take
  `thin`.** A model with `thin = t` runs its chain for `burnin + iterations * t`
  draws and keeps the last of every `t` after the burn-in. So `iterations` is
  still the number of draws kept and the posterior is still sized by it, while
  a slowly mixing chain can run long without every draw being held in memory. A
  chain thinned this way is exactly every `t`-th draw of the unthinned chain
  from the same seed, and a test asserts that for two dynamic factor models and
  the FAVAR. The `mcpar` of the draws and of the forecast says which were kept:
  iterations `t`, `2t`, ... after the burn-in. `model$thin` is carried only when
  it is above 1, so a model created without it is unchanged.

    **Draws are unchanged** for every model that does not thin: each binding
    reads a `thin` of 1, and the draws come back as the same `mcmc` objects as
    before.

* **Vendored BayesTS core refreshed: the samplers can thin.** **Draws are
  unchanged**, for all four dynamic factor models and the FAVAR, verified the
  same way as the refresh below. Each of the five was fitted from a pinned seed
  against the package built before and after it, through the coefficients, the
  log likelihood and a four-step forecast. Every posterior block of every one of
  them is bit-identical.

    Upstream's `VarSpec` gained `thin`, with `keeps()` and `kept_index()`, and
    every sampler now keeps its draws through those two calls instead of testing
    `draw >= burnin`. At `thin = 1` they are exactly that test. No binding here
    sets `thin` -- each fills `iterations` and `burnin` only -- so every model
    still keeps every draw after the burn-in. Exposing a thinning interval is a
    change for another time.

    The vendored set is still the same 32 files, so `inst/COPYRIGHTS` is
    unchanged.

* **Vendored BayesTS core refreshed again.** **Draws are unchanged**, for all
  four dynamic factor models and the FAVAR -- verified here rather than taken on
  trust. Each of the five was fitted from a pinned seed against the package built
  before and after the refresh, through the coefficients, the log likelihood and
  a four-step forecast, and every posterior block of every one of them is
  bit-identical.

    Nothing this package calls moved. Upstream made two changes, and both belong
    to the VAR and VEC samplers `bvartools` carries:
    - The out-of-sample regressors are now the compact layout, one period per
      row. `ForecastData::z` became `ForecastData::x`, and `bayests/spec.h`
      gained `VarSpec::n_x()`, which `n_non_structural()` is now written in
      terms of, to the same value. The two helpers in `core/models/model_support.h`
      that read those regressors changed with it. No factor sampler calls
      either one, and no binding here fills `ForecastData`.
    - A time varying cointegration space can be centred on a given space.
      `TvpCointSpacePrior` gained a `p_tau` transition, and `core/inputs.cpp`
      its validation. That is reached only from the error correction
      validators, which no factor model calls.

    The vendored set is still the same 32 files, so `inst/COPYRIGHTS` is
    unchanged.

* **Vendored BayesTS core refreshed.** **Draws are unchanged**, for all four
  dynamic factor models -- verified here rather than taken on trust: the four
  samplers were fingerprinted from a pinned seed before and after the refresh,
  and every posterior block of every one of them is bit-identical.

    Nothing this package calls moved. What upstream added is a prior on `rho`,
    the autocorrelation of a time varying cointegration space, which makes it a
    drawn parameter in the three time-varying error correction models
    `bvartools` carries. None of those is vendored here, and neither is the
    truncated normal the new block draws from. Three shared files change and
    all three are type surface: `bayests/priors.h` gains a `CointRhoPrior` and
    a field on `TvpCointSpacePrior`, `bayests/results.h` a `rho` member on the
    three `VecTvp*Draws`, and `core/inputs.cpp` the validation of the new
    prior's support -- reached only from the error correction validators, which
    no factor model calls.

    Upstream also corrected the state a time varying cointegration space is
    centred on in the first period, which moves the draws of those same three
    error correction models. Again none of them is vendored here.

* **`add_initial_values()` on class `dfmodel` starts the loadings at the
  principal components estimate, which keeps the sampler out of the mirror
  mode.** The loadings are the one block whose starting value decides which mode
  the chain converges to rather than how long it takes to get there. Only the
  product of the loadings and the factors is identified, and the restriction
  that pins it fixes the leading `n x n` block of lambda rather than anything
  about the factors, so a start whose free loadings are negative where the data
  want them positive is a coherent model in its own right -- the mirror, in
  which the factor is the negative of the common component and the series whose
  loading is fixed at one is treated as noise, its idiosyncratic variance
  absorbing nearly all of its variation. It fits far worse, and burn-in does not
  escape it, since turning the factor around would have to pass through
  configurations no single Gibbs step will take.

  The starting values were drawn from a prior whose default precision of 0.01 is
  a standard deviation of ten, so their signs were close to a coin toss and, on
  a ten-series panel of US real activity, roughly half of all seeds ended in the
  mirror -- on the constant model as readily as on the time-varying and
  stochastic volatility ones. The new default, `method = "pca"`, starts them at
  the loadings the leading `n` principal components of `x` imply, rotated so
  that the identifying block is the identity. That estimate is what the data
  say, it does not depend on the seed, and it is invariant to the sign
  convention `svd()` returns. Every other block is still drawn from its prior.

  `method = "prior"` keeps the old behaviour for whoever wants every block
  drawn, and an unsupported `method` is now an error rather than a model that
  silently comes back without starting values. The cause is what this removes
  rather than the possibility: on a sample too short to hold the chain, the
  mirror is still reachable, so the sign of the estimated loadings stays worth a
  look.

* **A vignette covers `error = "sv"` and `tvp = TRUE`.** *Drifting loadings and
  changing volatility* estimates all four dynamic factor model algorithms on one
  panel of US real activity: what the six-element stochastic volatility prior and
  the state equation on the coefficients ask for, how to read a posterior that
  holds one value per period rather than one value, what the two specifications
  do to each other when both are carried, and what they do to a forecast
  interval. Both vignettes are now precompiled from `.Rmd.orig` sources, as in
  `bvartools` and `bgvars`, since posterior simulation is too slow to run during
  `R CMD build`.

* **A dynamic factor model with a single observed series and a single factor is
  estimated rather than crashing.** It is the one specification whose loading
  matrix is entirely fixed -- the identifying restriction makes the leading
  `n x n` block of lambda unit lower triangular, which at `m = 1` and
  `n = 1` is the whole of it -- so the count of freely estimated loadings,
  `(2m - n - 1)n/2`, is zero. `add_priors()` built an empty prior for that block
  and `add_initial_values()` then called `chol()` on it, which reported a
  zero-dimensional matrix and said nothing about the model. Both now carry the
  empty block through, for constant and for time-varying coefficients alike, and
  all four samplers estimate the result: `x_t = f_t + u_t` with an AR(`p`)
  factor, which is an unobserved-components model.

* **`add_initial_values()` on class `dfmodel` draws its starting values at the
  scale of the prior rather than its inverse.** Every `vinv` here is a
  precision, so a draw from the prior has standard deviation `1/sqrt(vinv)`,
  which is what `backsolve(chol(vinv), z)` gives; the loadings and the
  transition coefficients of a constant-coefficient model instead used
  `chol(vinv) %*% z`, which gives `sqrt(vinv)`. The drifting blocks and the
  log-volatilities already had it right, so the file held both arrangements at
  once; all four now go through one function.

  The two forms agree at `vinv = 1` and nowhere else, and they diverge in
  opposite directions as the prior tightens: at the default precision of 0.01
  the starting values were a hundred times tighter than the prior, and at a
  precision of 100 they would have been a hundred times looser, so it was the
  caller stating a firm prior who was thrown furthest from it. This affects
  where the chains start and not what they converge to -- the sampler redraws
  every block from the data in its first iteration -- so it is burn-in rather
  than correctness.

* **`create_dfmodel()` accepts a univariate `ts`.** The branch meant to handle
  input with no column names assigned the matrix conversion to a local that
  nothing read, so `x` stayed a vector and naming it failed with `'dimnames'
  applied to non-array`. A `ts` matrix whose `dimnames` had been stripped reached
  the same branch and failed on the length of the names instead. Both are now
  converted and named -- `y` for a single column, `y1`, `y2`, ... for several --
  and the series keeps its time index, which `scale()` would otherwise have
  dropped from a vector. The condition keys on the column names rather than on
  `dimnames` as a whole, so a matrix carrying row names but no column names is
  named here too instead of carrying empty names onward.

* **`create_dfmodel()` rejects more factors than observed series**, naming the
  cause, instead of leaving it to the sampler. `n` equal to the number of columns
  of `x` is a valid specification and still runs -- the elements below the
  diagonal of the leading block are freely estimated -- but `n` above it is not a
  model, and the count of free loadings it implied was merely a wrong number that
  turned negative far enough out, where a caller saw `diag()` refusing a negative
  dimension. A vector `n` is checked entry by entry.

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
  a loading at zero to seven decimal places.

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
of the sampler. `bvartools (>= 0.3.0.9000)` is, for the generics
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
