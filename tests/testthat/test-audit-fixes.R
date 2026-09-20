# Three defects found in an audit of the package, and the shapes that pin them
# down. Each of these failed before the fix it is named for.

test_that("the error precisions start at a draw from their gamma prior", {

  # add_priors() stores a Gamma(shape, rate) prior on each precision and the
  # sampler draws each precision from that gamma updated by the data, so a
  # starting value is a draw from the gamma itself. Taking the reciprocal of a
  # gamma drawn with the rate inverted -- which this did -- puts the default
  # model twenty times below the prior mean of 5/4.
  sim <- sim_dfm(tt = 60, m = 6, n = 2)

  draws <- vapply(seq_len(400), function(i) {
    set.seed(1000 + i)
    object <- add_initial_values(add_priors(
      create_dfmodel(x = sim$x, p = 1, n = 2, iterations = 10, burnin = 5)))
    c(diag(object$initial$uinv), diag(object$initial$vinv))
  }, numeric(6 + 2))

  # Prior mean shape/rate = 5/4, standard deviation sqrt(shape)/rate = sqrt(5)/4.
  # 400 draws of each of eight elements put the standard error of a mean at
  # about 0.028, so three of them is a wide band that the reciprocal form
  # (mean 1/16) misses by a factor of twenty.
  expect_equal(mean(draws), 5 / 4, tolerance = 0.05)
  expect_gt(min(rowMeans(draws)), 1)
  expect_lt(max(rowMeans(draws)), 1.5)

  # A factor augmented VAR draws its idiosyncratic precisions from the same
  # helper, so the two cannot drift apart again.
  favar <- make_favar_sample(tt = 60, n_x = 5)
  favar_draws <- vapply(seq_len(400), function(i) {
    set.seed(2000 + i)
    object <- add_initial_values(add_priors(
      create_favarmodel(x = favar$x, y = favar$yts, p = 1, n = 1,
                        normalize_x = FALSE, iterations = 10, burnin = 5)))
    diag(object$initial$uinv)
  }, numeric(5))
  expect_equal(mean(favar_draws), 5 / 4, tolerance = 0.05)
})

test_that("an improper gamma prior on a precision is refused where it is written", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 10, burnin = 5)

  # rgamma() returns zero for a shape of zero, so nothing can be started from
  # this prior. Said by add_priors(), where the argument still has a name,
  # rather than one function later.
  expect_error(add_priors(object, u = list(shape = 0, rate = 4)),
               "'u[$]shape' must be larger than 0")
  expect_error(add_priors(object, v = list(shape = 0, rate = 4)),
               "'v[$]shape' must be larger than 0")

  # add_initial_values() keeps its own guard for a prior written by hand.
  hand_edited <- add_priors(object)
  hand_edited$priors$u$shape[] <- 0
  expect_error(add_initial_values(hand_edited), "must be larger than 0")
})

test_that("a cumulative response keeps its shape at horizon zero", {

  # apply() drops the horizon dimension where there is one of them, so the
  # transpose used to return a 1 x draws matrix and the quantiles were then
  # taken over the draws.
  est <- make_favar(iterations = 60, burnin = 40)
  object <- est$model

  plain <- irf(object, impulse = 1, response = 2, n_ahead = 0)
  cumulated <- irf(object, impulse = 1, response = 2, n_ahead = 0,
                   cumulative = TRUE)

  # One row, the impact horizon, and the three quantile columns.
  expect_equal(dim(cumulated), c(1L, 3L))
  # Accumulating a single horizon is the horizon itself.
  expect_equal(as.matrix(cumulated), as.matrix(plain), tolerance = 1e-12)

  draws <- as.matrix(irf(object, impulse = 1, response = 2, n_ahead = 0,
                         cumulative = TRUE, keep_draws = TRUE))
  expect_equal(ncol(draws), 1L)
  expect_equal(nrow(draws), nrow(as.matrix(
    irf(object, impulse = 1, response = 2, n_ahead = 0, keep_draws = TRUE))))
})

test_that("a scored test sample is put on the scale the model was estimated on", {

  # The forecast lives on the standardised panel, so realised values in the
  # data's own units have to be standardised with the estimation sample's
  # moments before they are scored against it. Scoring the raw values is
  # scoring the wrong quantity, and nothing about the result says so.
  sim <- sim_dfm(tt = 44)
  train <- stats::window(sim$x, end = c(1989, 4))
  test <- stats::window(sim$x, start = c(1990, 1))

  fit <- function(normalize_x) {
    object <- create_dfmodel(x = train, p = sim$p, n = sim$n,
                             normalize_x = normalize_x,
                             iterations = 20, burnin = 10)
    object <- add_initial_values(add_priors(object))
    set.seed(24)
    object <- add_posterior_coefficients(object)
    add_posterior_forecasts(object, n_ahead = 3)
  }

  scored <- add_predictive_loglik(fit(TRUE), test_sample = test)

  center <- attr(scored$data$x, "scaled:center")
  scale <- attr(scored$data$x, "scaled:scale")
  expected <- sweep(sweep(as.matrix(test)[1:3, , drop = FALSE], 2, center, "-"),
                    2, scale, "/")
  expect_equal(unname(scored$data$test$x), unname(expected))

  # Scoring the same model again, from what it carries, must not standardise a
  # second time.
  twice <- add_predictive_loglik(scored)
  expect_equal(twice$data$test$x, scored$data$test$x)
  expect_equal(unname(as.matrix(twice$posterior$forecast$loglik)),
               unname(as.matrix(scored$posterior$forecast$loglik)))

  # A model that did not normalise is scored in the data's units, untouched.
  raw <- add_predictive_loglik(fit(FALSE), test_sample = test)
  expect_equal(unname(raw$data$test$x), unname(as.matrix(test)[1:3, , drop = FALSE]))
})

test_that("the pointwise log-likelihood carries the labels of the kept draws", {

  # Every other block of draws is an mcmc object; this one was a bare matrix, so
  # a thinned model's log-likelihood said nothing about which draws it held.
  object <- prepared_dfm(iterations = 10, burnin = 5)$object
  object$model$thin <- 3L
  set.seed(5)
  object <- add_posterior_loglik(add_posterior_coefficients(object))

  loglik <- object$posterior$loglik
  expect_s3_class(loglik, "mcmc")
  expect_equal(attr(loglik, "mcpar"), c(3, 30, 3))
  expect_equal(nrow(loglik), 10L)

  favar <- make_favar_sample(tt = 60, n_x = 5)
  model <- create_favarmodel(x = favar$x, y = favar$yts, p = 1, n = 1,
                             normalize_x = FALSE, iterations = 10, burnin = 5,
                             thin = 2)
  model <- add_initial_values(add_priors(model))
  model <- add_posterior_loglik(add_posterior_coefficients(model))
  expect_s3_class(model$posterior$loglik, "mcmc")
  expect_equal(attr(model$posterior$loglik, "mcpar"), c(2, 20, 2))
})

test_that("create_dfmodel and create_favarmodel check their sampler arguments", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  favar <- make_favar_sample(tt = 60, n_x = 5)

  # These used to travel as far as the sampler, where the complaint arrived
  # through the error = TRUE path rather than as a stop, or -- for a character
  # p, which "a" < 0 compares as a string and lets through -- as far as a
  # matrix multiplication.
  expect_error(create_dfmodel(x = sim$x, p = "a", n = 1), "whole numbers")
  expect_error(create_dfmodel(x = sim$x, p = 1.5, n = 1), "whole numbers")
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 0),
               "at least 1")
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 2.7),
               "whole number")
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1, burnin = -3),
               "at least 0")

  expect_error(create_favarmodel(x = favar$x, y = favar$yts, p = 1, n = 1,
                                 iterations = 0), "at least 1")
  expect_error(create_favarmodel(x = favar$x, y = favar$yts, p = "a", n = 1),
               "whole numbers")

  # And the message a caller already got stays the one they get.
  expect_error(create_dfmodel(x = sim$x, p = -1, n = 1), "'p' must be at least 0")
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 0), "'n' must be at least 1")
})

test_that("add_initial_values says when the priors are missing", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 10, burnin = 5)
  expect_error(add_initial_values(object), "Did you call add_priors")
})
