# Scoring a forecast against what its horizon realised.
#
# A factor model is scored by filtering: each column of
# posterior$forecast$loglik is the density of that period's realised
# observation given the ones before it, the latent factors having been updated
# by each in turn. What the tests here pin is the shape it comes back in, that
# the realised values reach it, and that the model carries what it was scored
# against so the same call can be made again without them.

# A fitted model with a forecast, and the periods held back from it.
#
# `error` as well as `tvp`, so that all four samplers can be scored. Only the
# two constant-variance ones used to be: the filter behind the stochastic
# volatility score carries a log-volatility random walk over the horizon that
# nothing here ran, and the one bug this code has had was in that filter.
scored_setup <- function(n_ahead = 3, tvp = FALSE, error = "gamma",
                         iterations = 20, burnin = 10) {
  sim <- sim_dfm(tt = 44)
  train <- stats::window(sim$x, end = c(1989, 4))
  test <- stats::window(sim$x, start = c(1990, 1))

  object <- create_dfmodel(x = train, p = sim$p, n = sim$n, tvp = tvp,
                           error = error, iterations = iterations, burnin = burnin)

  coefficients <- if (tvp) {
    list(lambda = tvp_prior(vinv = 0.02, shape = 4, rate = 0.05),
         a = tvp_prior(vinv = 0.03, shape = 5, rate = 0.06))
  } else {
    list()
  }
  errors <- if (error == "sv") list(u = sv_prior(), v = sv_prior()) else list()
  object <- do.call(add_priors, c(list(object), coefficients, errors))

  object <- add_initial_values(object)
  set.seed(24)
  object <- add_posterior_coefficients(object)
  object <- add_posterior_forecasts(add_forecast_input(object, n_ahead = n_ahead))

  list(object = object, test = test, sim = sim)
}

test_that("every one of the four samplers can be scored", {

  # The two stochastic volatility scores were never run here. They work; what
  # was missing was anything that said so, which is how the defect in this
  # filter reached a release and was found by a toolchain elsewhere.
  for (error in c("gamma", "sv")) {
    for (tvp in c(FALSE, TRUE)) {
      fit <- scored_setup(tvp = tvp, error = error)
      scored <- add_predictive_loglik(fit$object, test_sample = fit$test)
      loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

      label <- paste0("error = ", error, ", tvp = ", tvp)
      expect_s3_class(loglik, "mcmc")
      expect_identical(dim(loglik), c(20L, 3L), label = label)
      expect_true(all(is.finite(loglik)), label = label)
    }
  }
})

test_that("a test sample the horizon never reaches leaves the model alone", {

  # Documented: "a model estimated to the end of a series forecasts past what
  # was ever observed", and that is not an error. The object comes back as it
  # went in, with no loglik and nothing in data$test.
  fit <- scored_setup()

  before <- stats::window(fit$sim$x, end = c(1980, 4))
  untouched <- add_predictive_loglik(fit$object, test_sample = before)

  expect_null(untouched[["posterior"]][["forecast"]][["loglik"]])
  expect_null(untouched[["data"]][["test"]])
  expect_identical(untouched, fit$object)

  # An empty sample is the same answer.
  empty <- stats::window(fit$sim$x, start = c(1980, 1), end = c(1980, 1))
  empty[] <- NA_real_
  expect_identical(add_predictive_loglik(fit$object, test_sample = empty), fit$object)
})

test_that("a forecast is scored against the periods its horizon realised", {
  fit <- scored_setup()
  scored <- add_predictive_loglik(fit$object, test_sample = fit$test)
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

  expect_s3_class(loglik, "mcmc")
  # One row per draw and one column per scored period, as every other block of
  # the posterior is laid out.
  expect_identical(dim(loglik), c(20L, 3L))
  expect_true(all(is.finite(loglik)))

  # The paths are still there beside it: the two are members of one group.
  expect_true(all(c("forecasts", "loglik") %in%
                    names(scored[["posterior"]][["forecast"]])))
})

test_that("a scored model carries what it was scored against", {
  fit <- scored_setup()
  scored <- add_predictive_loglik(fit$object, test_sample = fit$test)
  realised <- scored[["data"]][["test"]][["x"]]

  expect_identical(dim(realised), c(3L, as.integer(fit$sim$m)))

  # On the model's scale, not the data's: scored_setup() creates the model with
  # the default normalisation, so the forecast the density is computed against
  # is standardised and the realised values are put on that scale with the
  # estimation sample's own centre and spread before they are scored.
  train <- scored[["data"]][["x"]]
  expected <- stats::window(fit$test, end = c(1990, 3))
  expected <- sweep(sweep(as.matrix(expected), 2, attr(train, "scaled:center"), "-"),
                    2, attr(train, "scaled:scale"), "/")
  expect_equal(as.numeric(realised), as.numeric(expected))

  # And is scored again from them, to the same numbers, without the sample.
  again <- add_predictive_loglik(scored)
  expect_s3_class(again[["posterior"]][["forecast"]][["loglik"]], "mcmc")
  expect_equal(as.numeric(again[["posterior"]][["forecast"]][["loglik"]]),
               as.numeric(scored[["posterior"]][["forecast"]][["loglik"]]))
})

test_that("the realised values reach the density", {
  fit <- scored_setup()
  first <- add_predictive_loglik(fit$object, test_sample = fit$test)

  moved <- fit$test
  moved[1, ] <- moved[1, ] + 5
  second <- add_predictive_loglik(fit$object, test_sample = moved)

  a <- unclass(first[["posterior"]][["forecast"]][["loglik"]])
  b <- unclass(second[["posterior"]][["forecast"]][["loglik"]])

  expect_false(isTRUE(all.equal(a[, 1], b[, 1])))
  # And the column after it, which conditions on the period that moved: that is
  # what filtering means and what a score without the update would miss.
  expect_false(isTRUE(all.equal(a[, 2], b[, 2])))
})

test_that("fewer realised periods than the horizon are scored on their own", {
  fit <- scored_setup(n_ahead = 4)
  scored <- add_predictive_loglik(fit$object,
                                  test_sample = stats::window(fit$test, end = c(1990, 2)))

  expect_identical(ncol(scored[["posterior"]][["forecast"]][["loglik"]]), 2L)
})

test_that("a model with nothing to score against is refused", {
  fit <- scored_setup()
  expect_error(add_predictive_loglik(fit$object), "carries none in data")

  wide <- fit$object
  wide[["data"]][["test"]][["x"]] <- cbind(fit$test, fit$test[, 1])[1:3, , drop = FALSE]
  expect_error(add_predictive_loglik(wide), "columns, but the model has")

  long <- fit$object
  long[["data"]][["test"]][["x"]] <- as.matrix(fit$test)
  expect_error(add_predictive_loglik(long), "more than the 3 this model forecasts")
})

test_that("a model without a forecast cannot be scored", {
  sim <- sim_dfm(tt = 44)
  object <- create_dfmodel(x = sim$x, p = sim$p, n = sim$n, iterations = 20, burnin = 10)
  object <- add_posterior_coefficients(add_initial_values(add_priors(object)))

  expect_error(add_predictive_loglik(object, test_sample = sim$x), "no forecast horizon")
})

test_that("a model whose states drift is scored under them", {
  fit <- scored_setup(tvp = TRUE)
  scored <- add_predictive_loglik(fit$object, test_sample = fit$test)
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

  expect_identical(dim(loglik), c(20L, 3L))
  expect_true(all(is.finite(loglik)))

  # Holding the states is a different model from letting them drift, so it is a
  # different score. The states are drawn, hence the seed on either side.
  held <- fit$object
  held[["model"]][["forecast_states"]] <- "hold"
  set.seed(5)
  drifting <- add_predictive_loglik(fit$object, test_sample = fit$test)
  set.seed(5)
  holding <- add_predictive_loglik(held, test_sample = fit$test)

  expect_false(isTRUE(all.equal(
    unclass(drifting[["posterior"]][["forecast"]][["loglik"]]),
    unclass(holding[["posterior"]][["forecast"]][["loglik"]]))))
})
