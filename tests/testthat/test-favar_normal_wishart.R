test_that("the model object records the state's two halves separately", {
  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 2, n = 1,
                             iterations = 10, burnin = 5)

  expect_s3_class(model, "favarmodel")
  expect_equal(model$model$m, 8)
  expect_equal(model$model$n, 1)
  expect_equal(model$model$n_obs, 1)
  expect_equal(model$model$algorithm, "FavarNormalWishart")
  # The observed block is data of its own, not a column of the panel.
  expect_equal(ncol(model$data$y), 1)
  expect_equal(ncol(model$data$x), 8)
})

test_that("the priors have the widths the identification implies", {
  sim <- make_favar_sample()
  model <- add_priors(create_favarmodel(x = sim$x, y = sim$yts, p = 2, n = 1,
                                        iterations = 10, burnin = 5))

  n_state <- 2          # one factor plus one observed
  n_lambda <- (8 - 1) * n_state   # the first series carries no free loading

  expect_equal(dim(model$priors$lambda$vinv), c(n_lambda, n_lambda))
  expect_equal(dim(model$priors$a$vinv), c(n_state^2 * 2, n_state^2 * 2))
  expect_equal(nrow(model$priors$u$shape), 8)

  # Q is unrestricted, so its prior is a Wishart rather than independent gammas.
  expect_equal(model$priors$v$df, n_state)
  expect_equal(dim(model$priors$v$scale), c(n_state, n_state))
})

test_that("a Wishart prior too weak to be proper is refused", {
  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 2, n = 1,
                             iterations = 10, burnin = 5)
  expect_error(add_priors(model, v = list(df = 1, scale = NULL)),
               "at least the width of the state")
  expect_error(add_priors(model, v = list(df = NULL, scale = diag(1, 3))),
               "must be a 2 x 2 matrix")
})

test_that("the initial loadings are the rectangle the binding expects", {
  sim <- make_favar_sample()
  model <- add_initial_values(add_priors(
    create_favarmodel(x = sim$x, y = sim$yts, p = 2, n = 1,
                      iterations = 10, burnin = 5)))

  # One row per series that carries a loading, one column per state element.
  expect_equal(dim(model$initial$lambda), c(7, 2))
  expect_equal(dim(model$initial$vinv), c(2, 2))
  expect_equal(dim(model$initial$uinv), c(8, 8))
  expect_true(all(diag(model$initial$uinv) > 0))
})

test_that("more factors than panel series is refused", {
  sim <- make_favar_sample()
  expect_error(create_favarmodel(x = sim$x, y = sim$yts, n = 8),
               "must be smaller than the number of columns")
})

test_that("it recovers the factor and the identification survives", {
  fit <- make_favar()
  model <- fit$model
  expect_false(isTRUE(model$error))

  factors <- as.matrix(model$posterior$factors$coeffs)
  expect_equal(nrow(factors), 600)
  expect_equal(ncol(factors), length(fit$sim$f))

  # The sign of a factor is a convention here only up to the identification,
  # which pins it: series 1 loads +1, so the factor cannot come back flipped.
  expect_gt(stats::cor(colMeans(factors), fit$sim$f), 0.95)

  # The identification pins the scale as well as the sign: series 1 loads
  # exactly one, so the factor comes back on the units of the data rather than
  # on an arbitrary normalisation of its own.
  expect_lt(abs(stats::sd(colMeans(factors)) / stats::sd(fit$sim$f) - 1), 0.1)

  lambda <- matrix(colMeans(as.matrix(model$posterior$lambda$coeffs)), nrow = 8)
  expect_equal(dim(lambda), c(8, 2))

  # The identifying row: exactly one on the factor, exactly zero on the observed
  # block, in every draw rather than on average.
  all_draws <- as.matrix(model$posterior$lambda$coeffs)
  expect_true(all(all_draws[, 1] == 1))   # element (1, 1) of each k x n_state
  expect_true(all(all_draws[, 9] == 0))   # element (1, 2), column-major
})

test_that("the free loadings are found in the right places", {
  fit <- make_favar()
  lambda <- matrix(colMeans(as.matrix(fit$model$posterior$lambda$coeffs)), nrow = 8)

  # Loadings on the factor, over the rows that carry one.
  expect_gt(stats::cor(lambda[-1, 1], fit$sim$load_f[-1]), 0.95)
  expect_lt(max(abs(lambda[-1, 1] - fit$sim$load_f[-1])), 0.1)

  # And on the observed block, which is the half a dynamic factor model has no
  # place for. That these come back at all is the whole point of the model:
  # a DFM would push them into the idiosyncratic variance.
  expect_gt(stats::cor(lambda[-1, 2], fit$sim$load_y[-1]), 0.95)
  expect_lt(max(abs(lambda[-1, 2] - fit$sim$load_y[-1])), 0.1)
})

test_that("the forecast covers the observed block as well as the panel", {
  fit <- make_favar()
  model <- add_posterior_forecasts(fit$model, n_ahead = 6)

  expect_false(isTRUE(model$error))
  # Wider than the panel: each horizon carries k panel series followed by the
  # n_obs observed ones. This is the one forecast in the package that is not
  # h * k, and getting it wrong would silently truncate the observed block.
  expect_equal(ncol(model$posterior$forecast), 6 * (8 + 1))
  expect_equal(nrow(model$posterior$forecast), 600)
  expect_true(all(is.finite(model$posterior$forecast)))
})

test_that("the log likelihood is draws by periods", {
  fit <- make_favar()
  model <- add_posterior_loglik(fit$model)

  expect_false(isTRUE(model$error))
  expect_equal(dim(model$posterior$loglik), c(600, length(fit$sim$f)))
  expect_true(all(is.finite(model$posterior$loglik)))
})

test_that("a forecast without draws is refused", {
  sim <- make_favar_sample()
  model <- add_initial_values(add_priors(
    create_favarmodel(x = sim$x, y = sim$yts, p = 2, n = 1,
                      iterations = 10, burnin = 5)))
  expect_error(add_posterior_forecasts(model), "does not contain posterior draws")
})
