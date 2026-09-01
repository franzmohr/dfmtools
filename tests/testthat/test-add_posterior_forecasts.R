test_that("add_posterior_forecasts returns h x M columns per draw", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 2, p = 2)

  set.seed(21)
  object <- add_posterior_forecasts(add_posterior_coefficients(prep$object),
                                    n_ahead = 3)

  expect_s3_class(object, "dfmodel")
  expect_s3_class(object$posterior$forecast, "mcmc")
  expect_equal(dim(object$posterior$forecast), c(20L, 3L * 4L))
  expect_true(all(is.finite(object$posterior$forecast)))

  # The horizon is the whole of what the sampler needs: a dynamic factor model
  # has no out-of-sample regressors to supply.
  expect_equal(object$model$h, 3L)
})

test_that("the forecast columns are the variables within a horizon", {

  # One variable inflated a hundredfold and no normalisation, so the layout is
  # readable off the magnitudes: the documented order is the horizons stacked in
  # the variable order of the sample, i.e. columns m, 2m, ... belong to the last
  # variable.
  sim <- sim_dfm(tt = 40, m = 3, n = 1, p = 1)
  x <- sim$x
  x[, 3] <- x[, 3] * 100

  # Enough burn-in for the loadings to have found that scale; at 10 draws the
  # magnitudes still reflect the starting values rather than the data.
  object <- create_dfmodel(x = x, p = 1, n = 1, iterations = 100, burnin = 200,
                           normalize_x = FALSE)
  object <- add_initial_values(add_priors(object))

  set.seed(25)
  object <- add_posterior_forecasts(add_posterior_coefficients(object), n_ahead = 4)

  scale <- matrix(colMeans(abs(object$posterior$forecast)), nrow = 3, ncol = 4)

  # Row 3 of that M x h reshape is the inflated variable, and it is the largest
  # of the three at every horizon. Note that variable 1 is not a useful
  # comparison: its loading is fixed at 1 by the identifying restriction, so it
  # follows whatever scale the factor takes on.
  expect_equal(apply(scale, 2, which.max), rep(3L, 4))
  expect_true(all(scale[3, ] > 10 * scale[2, ]))
})

test_that("forecasts are reproducible under set.seed", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)

  set.seed(22)
  drawn <- add_posterior_coefficients(prep$object)

  set.seed(23)
  first <- add_posterior_forecasts(drawn, n_ahead = 4)$posterior$forecast

  set.seed(23)
  second <- add_posterior_forecasts(drawn, n_ahead = 4)$posterior$forecast

  expect_equal(first, second)
})

test_that("add_posterior_forecasts rejects invalid input", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)

  expect_error(add_posterior_forecasts(prep$object, n_ahead = 4),
               "does not contain posterior draws")

  set.seed(24)
  drawn <- add_posterior_coefficients(prep$object)

  expect_error(add_posterior_forecasts(drawn, n_ahead = 0),
               "'n_ahead' must be at least 1")

  # The factor path cannot be recomputed from the parameters, so its absence is
  # an error rather than an extra filtering pass.
  without_factors <- drawn
  without_factors$posterior$factors <- NULL
  expect_error(add_posterior_forecasts(without_factors, n_ahead = 4),
               "does not contain posterior draws of the factors")
})
