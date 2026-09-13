---
name: dfmtools
description: Write correct dfmtools code — Bayesian dynamic factor models and factor augmented VARs in R with create_dfmodel, create_favarmodel, add_priors, add_initial_values, add_posterior_coefficients, add_posterior_forecasts, add_posterior_loglik, irf and fevd. Use whenever R code loads dfmtools, handles a 'dfmodel' or 'favarmodel' object, or estimates a DFM, a FAVAR, factor loadings, a factor path, or a factor model with time-varying loadings or stochastic volatility in R. Also use when reading loadings, factors or forecasts out of such an object, or when identifying a shock in a FAVAR.
---

# Writing correct dfmtools code

dfmtools estimates dynamic factor models (DFM)

    x_t = lambda f_t + u_t,    f_t = A_1 f_{t-1} + ... + A_p f_{t-p} + v_t

and factor augmented VARs (FAVAR), where observed variables `y_t` join the
factors in the state. It depends on `bvartools`: `library(dfmtools)` attaches
it, and the steps are bvartools generics with `dfmodel` and `favarmodel` methods.
The samplers are the C++ core of BayesTS, and `set.seed()` reaches them.

```r
library(dfmtools)
set.seed(42)

data("bem_dfmdata")                        # 225 quarters of 196 US series
x <- bem_dfmdata[, 1:20]

model <- create_dfmodel(x = x, p = 1, n = 2, iterations = 500, burnin = 200)
model <- add_priors(model,
                    lambda = list(vinv = 0.01), u = list(shape = 5, rate = 4),
                    a = list(vinv = 0.01), v = list(shape = 5, rate = 4))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)
model <- add_posterior_forecasts(model, n_ahead = 4)
```

`iterations = 500` keeps examples quick. The defaults are 20000 after 2000
burn-in, and a real analysis needs that order. `thin = t` runs `t` times as long
and keeps the last of every `t` draws, so a slowly mixing chain can run long
while the posterior still holds `iterations` draws. The draws that are not kept
are never stored, and `coda::mcpar()` counts the draws the chain actually ran:

```r
long <- create_dfmodel(x = x, p = 1, n = 2, iterations = 50, burnin = 50, thin = 5)
long <- add_priors(long,
                   lambda = list(vinv = 0.01), u = list(shape = 5, rate = 4),
                   a = list(vinv = 0.01), v = list(shape = 5, rate = 4))
long <- add_initial_values(long)
long <- add_posterior_coefficients(long)

stopifnot(nrow(long$posterior$lambda$coeffs) == 50,
          all(coda::mcpar(long$posterior$lambda$coeffs) == c(5, 250, 5)))
```

Reach for it when a posterior summary moves with the seed: the chain is too short
for how slowly it mixes.

## The rules that prevent wrong results

**1. Assign every step back**: `model <- add_priors(model, ...)`.

**2. The order of `x` is the identification.** Only `lambda f_t` is identified.
A DFM fixes the leading `n x n` block of `lambda` to be unit lower triangular,
so the first `n` columns of `x` define the factors, and `n` must not exceed the
number of series. Put broad, well-measured series first. A FAVAR fixes that
block to the **identity** instead, so factor `i` *is* panel series `i` up to
noise. The two rules are not interchangeable, and neither is the number of free
loadings they leave.

**3. Check the sign of the loadings.** A start whose free loadings have the
wrong sign is a coherent mirror model, and the sampler stays in it. The default
`add_initial_values(method = "pca")` starts from principal components to avoid
that, but on a short or weakly correlated panel the mirror is still reached.
Loadings that are negative where the panel is positively correlated are the
symptom.

**4. `x` is normalised by default.** `normalize_x = TRUE` scales every column,
so loadings, forecasts and panel impulse responses are in standard deviations of
each series. `attr(model$data$x, "scaled:scale")` and `"scaled:center"` undo it.

**5. Draws are rows.** Every `model$posterior$<block>$coeffs` is a `coda::mcmc`
matrix, draws × parameters.
- `lambda$coeffs` holds the **whole** `m x n` loading matrix per draw, fixed
  block included, column-major: `matrix(colMeans(lambda), nrow = m)`. For a
  FAVAR it is `k x (n + n_obs)`.
- `factors$coeffs` stores period by period, all `n` factors of a period
  together, so factor `i` is `f[, seq(i, tt * n, by = n)]`.
- `forecast` has `h × m` columns, horizon by horizon with the variables of each
  horizon in sample order. A FAVAR's covers `h × (k + n_obs)`, panel then
  observed.
- `loglik` is draws × periods, and it is **conditional on the drawn factor
  path**, not the marginal likelihood.

**6. Forecasts take the horizon directly.** `add_posterior_forecasts(model,
n_ahead = h)`. There are no out-of-sample regressors, and `add_forecast_input()`
has no method for these models.

**7. The prior names are not bvartools'.** Loadings and transition take
`vinv`, a precision: `lambda = list(vinv = 0.01)`, not `v_i`. A wrong name is an
error. Unlike bvartools, `add_priors()` has working defaults for the gamma
model. The stochastic volatility lists do use bvartools' names:
- `tvp = TRUE`: `lambda` and `a` also need `shape` and `rate` for the random walk.
- `error = "sv"`: `u` and `v` need `mu`, `v_i`, `shape`, `rate`,
  `state_variance` and `offset`.
- FAVAR: `v` is a Wishart prior (`df`, `scale`), because `Q` is unrestricted.

**8. There is no `irf()`, `fevd()` or `selection_criteria()` for a DFM.** Only a
FAVAR has impulse responses and variance decompositions. Model comparison on a
DFM works from `add_posterior_loglik()` directly.

**9. Vectors make lists.** `p = 1:2` or `n = 1:3` give a `'modellist'` that
every step maps over.

## FAVAR specifics

`create_favarmodel(x, y, p, n)` takes the panel `x` and the observed variables
`y`. **`y` is part of the state, not a set of regressors**: the transition is a
VAR over `(f_t, y_t)`, and `Q` is estimated in full so that its cross block can
be read.

- `add_priors(model, slow = c(...))` names panel series assumed not to respond
  to `y` within the period. Their loadings on `y` are held at zero, which is the
  Bernanke, Boivin and Eliasz restriction.
- `irf(model, impulse, response, n_ahead)`: `impulse` is an element of the
  **state** (a factor or an observed variable), and `response` is any observable,
  the panel then `y`. Factor `i` is named after panel column `i`.
- `type = "oir"` (default) uses the Cholesky factor of `Q` in `order`, and that
  ordering is an assumption the model cannot test. `"gir"` answers a different
  question and undoes the `slow` restriction on impact.

## Reading the installed documentation

The help pages are the reference. Read the method, not the generic —
`add_priors.dfmodel` and `add_priors.favarmodel` differ:

```r
txt <- capture.output(tools::Rd2txt(
  tools::Rd_db("dfmtools")[["add_priors.dfmodel.Rd"]],
  options = list(underline_titles = FALSE)))
stopifnot(length(txt) > 0)
```

Vignettes: `vignette("tvp-and-stochastic-volatility", package = "dfmtools")` and
`vignette("favar-monetary-policy", package = "dfmtools")`. The copy of this
skill matching the installed version is at
`system.file("agents", package = "dfmtools")`.

## Reference files

| File | Contents |
| --- | --- |
| `references/recipes.md` | Complete, tested examples: a DFM with loadings, factors and forecasts read back; a TVP-SV DFM; a FAVAR with an impulse response |

For VAR and VEC models, see the `bvartools` skill.
