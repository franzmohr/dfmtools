# Deterministic terms in the factor models, and the forecast route that reads
# them: add_forecast_input() -> add_posterior_forecasts() -> predict(), as in
# bvartools.

# A panel with a factor structure on top of a level and a seasonal pattern that
# differ from series to series, so that C has something to find. Quarterly,
# starting in the second quarter, so that the seasonal dummies cannot line up
# with the first row by accident.
sim_seasonal_dfm <- function(tt = 120, m = 5, seed = 7) {
  sim <- sim_dfm(tt = tt, m = m, seed = seed)
  season <- stats::cycle(stats::ts(seq_len(tt), start = c(1980, 2), frequency = 4))
  c_true <- cbind(const = seq(0.5, 2.5, length.out = m),
                  season.1 = rep(c(1.5, -1.5), length.out = m),
                  season.2 = 0,
                  season.3 = rep(c(-1, 1), length.out = m))
  d <- cbind(1, season == 1, season == 2, season == 3)
  x <- unclass(sim$x) + d %*% t(c_true)
  list(x = stats::ts(x, start = c(1980, 2), frequency = 4), c = c_true,
       lambda = sim$lambda)
}

test_that("the terms are bvartools' terms, named and valued alike", {
  x <- stats::ts(matrix(stats::rnorm(60), 20, 3, dimnames = list(NULL, c("a", "b", "c"))),
                 start = c(2000, 2), frequency = 4)

  model <- create_dfmodel(x, p = 1, n = 1, deterministic = "both", seasonal = TRUE,
                          iterations = 10, burnin = 10)
  var <- bvartools::create_bvarmodel(x, p = 0, deterministic = "both", seasonal = TRUE,
                                     iterations = 10, burnin = 10)

  expect_identical(model$model$deterministic, c("const", "trend", "season.1", "season.2", "season.3"))
  expect_equal(unclass(as.matrix(model$data$deterministic)),
               unclass(as.matrix(var$data$original$deterministic)), ignore_attr = TRUE)
  expect_identical(stats::tsp(model$data$deterministic), stats::tsp(x))
})

test_that("a model without them carries nothing of them", {
  prep <- prepared_dfm()
  expect_null(prep$object$model$deterministic)
  expect_null(prep$object$data$deterministic)
  expect_null(prep$object$priors$c)
  expect_null(prep$object$initial$c)
})

test_that("the specification is checked where it is made", {
  x <- sim_dfm()$x
  expect_error(create_dfmodel(x, p = 1, n = 1, deterministic = "level"), "deterministic")
  expect_error(create_dfmodel(x, p = 1, n = 1, seasonal = TRUE), "'const' or 'both'")
  expect_error(create_favarmodel(x[, -1], x[, 1, drop = FALSE], p = 1, n = 1,
                                 deterministic = "trend", seasonal = TRUE),
               "'const' or 'both'")
})

test_that("the horizon continues the sample's terms", {
  sim <- sim_seasonal_dfm(tt = 21)
  model <- create_dfmodel(sim$x, p = 1, n = 1, deterministic = "both", seasonal = TRUE,
                          iterations = 10, burnin = 10)
  fcst <- prepare_forecast_input(model, n_ahead = 6)

  # Built over the sample and the horizon at once, the horizon is the tail.
  span <- stats::ts(matrix(0, 27, 1), start = c(1980, 2), frequency = 4)
  whole <- create_dfmodel(stats::ts(matrix(stats::rnorm(27 * 5), 27, 5), start = c(1980, 2),
                                    frequency = 4),
                          p = 1, n = 1, deterministic = "both", seasonal = TRUE,
                          iterations = 10, burnin = 10)$data$deterministic
  expect_identical(fcst$h, 6L)
  expect_equal(fcst$x, unclass(as.matrix(whole))[22:27, ], ignore_attr = TRUE)
  expect_identical(colnames(fcst$x), model$model$deterministic)

  # And a model without them needs nothing but the horizon.
  expect_null(prepare_forecast_input(prepared_dfm()$object, n_ahead = 3)$x)
})

test_that("given terms must cover the horizon", {
  sim <- sim_seasonal_dfm(tt = 21)
  model <- create_dfmodel(sim$x, p = 1, n = 1, deterministic = "const",
                          iterations = 10, burnin = 10)
  ok <- stats::ts(matrix(2, 3, 1), start = c(1985, 3), frequency = 4)
  expect_equal(prepare_forecast_input(model, n_ahead = 3, deterministic = ok)$x,
               matrix(2, 3, 1, dimnames = list(NULL, "const")))
  late <- stats::ts(matrix(2, 3, 1), start = c(1985, 4), frequency = 4)
  expect_error(prepare_forecast_input(model, n_ahead = 3, deterministic = late), "must cover")
  expect_error(prepare_forecast_input(prepared_dfm()$object, n_ahead = 3, deterministic = ok),
               "no deterministic terms")
})

test_that("C is recovered and the forecast carries it", {
  sim <- sim_seasonal_dfm()
  set.seed(11)
  model <- create_dfmodel(sim$x, p = 1, n = 1, normalize_x = FALSE,
                          deterministic = "const", seasonal = TRUE,
                          iterations = 600, burnin = 600)
  model <- add_priors(model)
  model <- add_initial_values(model)

  # Least squares, under either method, so nothing about it is random.
  expect_length(model$initial$c, 5 * 4)
  expect_equal(nrow(model$priors$c$vinv), 5 * 4)

  model <- add_posterior_coefficients(model)
  expect_null(model$error)
  draws <- model$posterior$c$coeffs
  expect_s3_class(draws, "mcmc")
  expect_identical(colnames(draws)[c(1, 6)], c("x1.const", "x1.season.1"))

  c_mean <- matrix(colMeans(draws), 5, 4)
  # The seasonal columns trade off against nothing; the constant against the
  # level of a persistent factor, which is why it is held to less.
  expect_lt(max(abs(c_mean[, 2:4] - sim$c[, 2:4])), 0.4)
  expect_lt(max(abs(c_mean[, 1] - sim$c[, 1])), 1)

  model <- add_forecast_input(model, n_ahead = 8)
  model <- add_posterior_forecasts(model)
  pred <- predict(model)
  expect_s3_class(pred, "bvarprd")
  expect_identical(dim(pred$fcst), c(8L, 5L, 600L))

  # Far enough out, the mean of the forecast is C d_{T+h}, the seasonal pattern
  # among it.
  want <- model$data$forecast$x %*% t(c_mean)
  got <- apply(pred$fcst, 1:2, mean)
  expect_lt(max(abs(got[5:8, ] - want[5:8, ])), 0.4)
})

test_that("a forecast of a model with them needs their horizon", {
  sim <- sim_seasonal_dfm(tt = 40)
  model <- create_dfmodel(sim$x, p = 1, n = 1, deterministic = "const",
                          iterations = 20, burnin = 10)
  model <- add_posterior_coefficients(add_initial_values(add_priors(model)))
  model$model$h <- 4L
  expect_error(add_posterior_forecasts(model), "add_forecast_input")
})

test_that("every factor model estimates and forecasts with them", {
  sim <- sim_seasonal_dfm(tt = 40)
  specs <- list(list(error = "gamma", tvp = FALSE), list(error = "sv", tvp = FALSE),
                list(error = "gamma", tvp = TRUE), list(error = "sv", tvp = TRUE))
  for (spec in specs) {
    model <- create_dfmodel(sim$x, p = 1, n = 1, error = spec$error, tvp = spec$tvp,
                            deterministic = "both", seasonal = TRUE,
                            iterations = 20, burnin = 10)
    priors <- list(model)
    if (spec$tvp) {
      priors$lambda <- tvp_prior()
      priors$a <- tvp_prior()
    }
    if (spec$error == "sv") {
      priors$u <- sv_prior()
      priors$v <- sv_prior()
    }
    model <- add_initial_values(do.call(add_priors, priors))
    model <- add_posterior_coefficients(model)
    expect_null(model$error, label = model$model$algorithm)
    # A constant, a trend and three seasonal dummies for each of five series.
    expect_identical(dim(model$posterior$c$coeffs), c(20L, 5L * 5L))

    model <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 3))
    model <- add_posterior_loglik(model)
    model <- add_predictive_loglik(model, test_sample = sim$x[38:40, ] + 0)
    expect_true(all(is.finite(model$posterior$forecast$forecasts)), label = model$model$algorithm)
    expect_true(all(is.finite(model$posterior$loglik)), label = model$model$algorithm)
  }
})

test_that("the old way of setting the horizon still works, with a warning", {
  prep <- prepared_dfm()
  model <- add_posterior_coefficients(prep$object)

  set.seed(3)
  expect_warning(old <- add_posterior_forecasts(model, n_ahead = 3), "deprecated")
  set.seed(3)
  new <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 3))
  expect_identical(old$posterior$forecast$forecasts, new$posterior$forecast$forecasts)

  set.seed(3)
  expect_warning(default <- add_posterior_forecasts(model), "set to 10")
  expect_identical(default$model$h, 10L)
})

test_that("predict undoes the normalisation", {
  prep <- prepared_dfm()
  model <- add_posterior_coefficients(prep$object)
  model <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 2))

  pred <- predict(model)
  expect_equal(unclass(as.matrix(pred$y)), unclass(as.matrix(prep$sim$x)), ignore_attr = TRUE)
  expect_identical(stats::tsp(pred$y), stats::tsp(prep$sim$x))

  raw <- predict(model, rescale = FALSE)
  scale <- attr(model$data$x, "scaled:scale")
  center <- attr(model$data$x, "scaled:center")
  expect_equal(pred$fcst[, 2, ], raw$fcst[, 2, ] * scale[2] + center[2])
  expect_warning(predict(model, n_ahead = 5), "Limiting")
})

test_that("forecast errors are realised less forecast, on the model's scale", {
  prep <- prepared_dfm(tt = 44)
  sim <- prep$sim
  train <- stats::window(sim$x, end = c(1989, 4))
  model <- create_dfmodel(train, p = 1, n = 1, iterations = 20, burnin = 10)
  model <- add_posterior_coefficients(add_initial_values(add_priors(model)))
  model <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 6))

  model <- add_forecast_errors(model, test_sample = sim$x)
  errors <- get_forecast_errors(model)
  expect_s3_class(errors, "mcmc")
  # Four realised periods of the six forecast.
  expect_identical(dim(errors), c(20L, 4L * 4L))
  realised <- model$data$test$x
  expect_equal(unclass(errors)[1, 1:4],
               realised[1, ] - unclass(model$posterior$forecast$forecasts)[1, 1:4],
               ignore_attr = TRUE)

  # Scored again from what it carries.
  again <- add_forecast_errors(model)
  expect_identical(get_forecast_errors(again), errors)
})

test_that("a FAVAR's observed variables get a mean of their own", {
  sample <- make_favar_sample(tt = 160)
  y <- sample$yts + 3
  set.seed(5)
  model <- create_favarmodel(sample$x, y, p = 1, n = 1, deterministic = "const",
                             iterations = 400, burnin = 400)
  model <- add_initial_values(add_priors(model))
  expect_length(model$initial$c_obs, 1)
  model <- add_posterior_coefficients(model)
  expect_null(model$error)
  expect_identical(colnames(model$posterior$c_obs$coeffs), "obs.const")

  # The mean the observed variable deviates from is what it has.
  expect_lt(abs(mean(model$posterior$c_obs$coeffs) - mean(y)), 0.5)

  model <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 4))
  pred <- predict(model)
  names <- c(colnames(model$data$x), "obs")
  expect_identical(dimnames(pred$fcst)[[2]], names)

  # Errors against a sample holding the panel and the observed variable, put
  # into the forecast's order by name.
  test <- cbind(y, sample$x)
  colnames(test) <- c("obs", colnames(model$data$x))
  short_x <- stats::window(sample$x, end = c(2026, 2))
  short_y <- stats::window(y, end = c(2026, 2))
  shortened <- create_favarmodel(short_x, short_y, p = 1, n = 1, deterministic = "const",
                                 iterations = 20, burnin = 10)
  shortened <- add_posterior_coefficients(add_initial_values(add_priors(shortened)))
  shortened <- add_posterior_forecasts(add_forecast_input(shortened, n_ahead = 4))
  shortened <- add_forecast_errors(shortened, test_sample = test)
  expect_identical(dim(get_forecast_errors(shortened)), c(20L, 4L * 9L))
  expect_identical(colnames(shortened$data$test$y), "obs")
})
