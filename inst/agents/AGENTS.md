# AGENTS.md — using dfmtools

Guidance for coding agents writing R code with
[dfmtools](https://github.com/franzmohr/dfmtools), which estimates Bayesian
dynamic factor models (DFM) and factor augmented VARs (FAVAR). The full skill is
`skills/dfmtools/SKILL.md`, with tested examples in
`skills/dfmtools/references/recipes.md`. This file is the part worth keeping in
context at all times.

1. **Assign every step back**: `model <- add_priors(model, ...)`.
2. **The order**: `create_dfmodel()` or `create_favarmodel()`, then
   `add_priors()`, `add_initial_values()`, `add_posterior_coefficients()`,
   `add_posterior_forecasts(model, n_ahead = h)` and `add_posterior_loglik()`.
   There is no `add_forecast_input()` step.
3. **Column order of `x` identifies the factors.** A DFM fixes the leading
   `n x n` block of the loadings to be unit lower triangular. A FAVAR fixes it to
   the identity, so its first `n` panel series *are* the factors up to noise.
4. **Look at the sign of the loadings.** A mirror mode with the factor negated
   is a coherent model the sampler can get stuck in. Keep the default
   `add_initial_values(method = "pca")`.
5. **`x` is standardised by default** (`normalize_x = TRUE`). Loadings, forecasts
   and responses are in standard deviations. Undo it with
   `attr(model$data$x, "scaled:scale")` and `"scaled:center"`.
6. **Draws are rows**: `model$posterior$<block>$coeffs` is draws × parameters.
   The loadings are stored whole (`m × n` per draw, column-major), the factor path
   period by period, and forecasts horizon by horizon.
7. **Prior precisions are called `vinv`** for `lambda` and `a`, not bvartools'
   `v_i`. The stochastic volatility lists do use bvartools' names.
8. **In a FAVAR, `y` is state, not regressors.** `irf()` shocks a state element
   and reads any observable, and a Cholesky ordering is an untestable assumption.
   `add_priors(model, slow = ...)` imposes the slow-moving restriction.
9. **No `irf()`, `fevd()` or `selection_criteria()` for a DFM**, only for a
   FAVAR.
10. **Keep `set.seed()`**: it reaches the C++ samplers. If two seeds give
    different summaries, the chain is too short: create the model with `thin`
    to run it longer without storing more draws.

The installed package carries these files at
`system.file("agents", package = "dfmtools")`. VAR and VEC models are in
`bvartools`, which has its own guide.
