# Worked examples

Every `r` block on this page runs in the package's test suite
(`tests/testthat/test-agent-docs.R`), top to bottom in one session. The shapes
asserted with `stopifnot()` are therefore what the installed version produces.
The draw counts are kept small so the tests stay quick. A real analysis needs
thousands of draws.

## A dynamic factor model

```r
library(dfmtools)
set.seed(42)

data("bem_dfmdata")
x <- bem_dfmdata[, 1:20]           # the first n columns define the factors

model <- create_dfmodel(x = x, p = 1, n = 2, iterations = 500, burnin = 200)
model <- add_priors(model,
                    lambda = list(vinv = 0.01),        # precision of the free loadings
                    u = list(shape = 5, rate = 4),     # idiosyncratic error precisions
                    a = list(vinv = 0.01),             # precision of the transition
                    v = list(shape = 5, rate = 4))     # factor innovation precisions
model <- add_initial_values(model)                     # method = "pca"
model <- add_posterior_coefficients(model)

m <- model$model$m
n <- model$model$n
tt <- nrow(model$data$x)
stopifnot(m == 20, n == 2)
```

The loadings come back whole, the fixed block included:

```r
lambda <- model$posterior$lambda$coeffs
stopifnot(inherits(lambda, "mcmc"), all(dim(lambda) == c(500, m * n)))

L <- matrix(colMeans(lambda), nrow = m, ncol = n)
# The identifying block: unit lower triangular in the first n rows
stopifnot(isTRUE(all.equal(L[1:n, 1:n][upper.tri(diag(n), diag = TRUE)], c(1, 0, 1))))
```

The factor path is stored period by period:

```r
f <- model$posterior$factors$coeffs
stopifnot(all(dim(f) == c(500, tt * n)))

factor_1 <- f[, seq(1, tt * n, by = n)]    # draws x periods
stopifnot(ncol(factor_1) == tt)
```

Forecasts need only the horizon:

```r
model <- add_posterior_forecasts(model, n_ahead = 4)
fc <- model$posterior$forecast
stopifnot(all(dim(fc) == c(500, 4 * m)))

series_1 <- fc[, seq(1, 4 * m, by = m)]    # series 1 at horizons 1 to 4

# In standard deviations of each series; back to the data's units:
scale_1 <- attr(model$data$x, "scaled:scale")[1]
center_1 <- attr(model$data$x, "scaled:center")[1]
series_1_units <- series_1 * scale_1 + center_1

refused <- tryCatch(add_forecast_input(model, n_ahead = 4), error = function(e) "no method")
stopifnot(identical(refused, "no method"))
```

The log likelihood is draws × periods, conditional on the factor path:

```r
model <- add_posterior_loglik(model)
stopifnot(all(dim(model$posterior$loglik) == c(500, tt)))
```

## Drifting loadings and changing volatility

`tvp = TRUE` gives `lambda` and `a` a random walk each, with `shape` and `rate`.
`error = "sv"` gives `u` and `v` the six stochastic volatility elements:

```r
tvp_sv <- create_dfmodel(x = x, p = 1, n = 1, tvp = TRUE, error = "sv",
                         iterations = 300, burnin = 100)
tvp_sv <- add_priors(tvp_sv,
                     lambda = list(vinv = 0.01, shape = 3, rate = 0.01),
                     a = list(vinv = 0.01, shape = 3, rate = 0.01),
                     u = list(mu = 0, v_i = 0.1, shape = 3, rate = 0.2,
                              state_variance = 0.05, offset = 1e-4),
                     v = list(mu = 0, v_i = 0.1, shape = 3, rate = 0.2,
                              state_variance = 0.05, offset = 1e-4))
tvp_sv <- add_initial_values(tvp_sv)
tvp_sv <- add_posterior_coefficients(tvp_sv)

# The loadings are now a path: the whole m x 1 matrix in every period
stopifnot(ncol(tvp_sv$posterior$lambda$coeffs) == m * 1 * tt)
```

## A factor augmented VAR

```r
panel <- bem_dfmdata[, 2:21]
observed <- bem_dfmdata[, 1, drop = FALSE]    # enters the state, not the regressors

favar <- create_favarmodel(x = panel, y = observed, p = 1, n = 2,
                           iterations = 500, burnin = 200)
favar <- add_priors(favar)                    # v is a Wishart prior: Q is unrestricted
favar <- add_initial_values(favar)
favar <- add_posterior_coefficients(favar)

k <- ncol(panel)
n_obs <- ncol(observed)
n_state <- 2 + n_obs
stopifnot(all(dim(favar$posterior$lambda$coeffs) == c(500, k * n_state)))

favar <- add_posterior_forecasts(favar, n_ahead = 4)
stopifnot(all(dim(favar$posterior$forecast) == c(500, 4 * (k + n_obs))))
```

An impulse is an element of the state, and a response any observable. Under
`type = "oir"` the Cholesky ordering in `order` is an assumption:

```r
ir <- irf(favar, impulse = colnames(observed), response = colnames(observed),
          n_ahead = 8)
stopifnot(inherits(ir, "bvarirf"), nrow(ir) == 9)
```

A dynamic factor model has no impulse responses of its own:

```r
refused <- tryCatch(irf(model, impulse = 1, response = 1, n_ahead = 4),
                    error = function(e) "no method")
stopifnot(identical(refused, "no method"))
```
