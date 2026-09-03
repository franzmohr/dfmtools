# The time varying parameter model, end to end through the R interface.
#
# What the numerics do is covered upstream, in BayesTS's
# test/unit_dfm_tvp_gamma.cpp -- the stacking conventions, the factor path
# against a dense posterior, the recovery of a loading that moves. What is left
# for here is the translation, and this model adds three things to it that the
# constant-coefficient ones do not have:
#
#   * two coefficient blocks that are paths rather than points, each with a state
#     equation of its own. The prior group gains `shape` and `rate`, the initial
#     values gain a path, a state before the sample and a precision, and none of
#     that reaches the sampler unless it is named correctly on both sides.
#   * a posterior whose `lambda` and `a` are tt times wider than before, plus a
#     `sigma` element that did not exist. Anything downstream that slices them by
#     the old width is silently wrong rather than broken.
#   * two entry points that want different slices of it. The forecast holds the
#     coefficients at their last in-sample period and the log likelihood scores
#     every period under its own, so the binding cuts one and not the other.
#
# One recovery test at the end, because the translation being right is not the
# same as the model finding a loading that moves, and the R side is where a
# permuted prior or a transposed path would show up as a plausible number.

test_that("create_dfmodel selects the time varying sampler", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1, p = 1)

  object <- create_dfmodel(x = sim$x, p = 1, n = 1, tvp = TRUE,
                           iterations = 20, burnin = 10)

  expect_s3_class(object, "dfmodel")
  expect_true(object$model$tvp)
  expect_equal(object$model$error, "gamma")
  expect_equal(object$model$algorithm, "DfmTvpGamma")

  # And the constant-coefficient default is untouched, tvp included.
  constant <- create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 20, burnin = 10)
  expect_false(constant$model$tvp)
  expect_equal(constant$model$algorithm, "DfmNormalGamma")

  # Refused rather than silently ignored: there is no sampler for it.
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv", tvp = TRUE,
                              iterations = 20, burnin = 10),
               "only available for error")

  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1, tvp = "yes",
                              iterations = 20, burnin = 10),
               "must be of class 'logical'")
})

test_that("add_priors builds a state equation for both coefficient blocks", {

  sim <- sim_dfm(tt = 40, m = 5, n = 2, p = 2)
  object <- create_dfmodel(x = sim$x, p = 2, n = 2, tvp = TRUE,
                           iterations = 20, burnin = 10)

  object <- add_priors(object,
                       lambda = tvp_prior(vinv = 0.02, shape = 4, rate = 0.05),
                       a = tvp_prior(vinv = 0.03, shape = 5, rate = 0.06))

  n_lambda <- n_free_lambda(5, 2)
  n_a <- 2 * 2 * 2

  # `vinv` keeps its meaning as a precision matrix and gains the pair beside it.
  expect_equal(dim(object$priors$lambda$vinv), rep(n_lambda, 2))
  expect_equal(diag(object$priors$lambda$vinv), rep(0.02, n_lambda))
  expect_equal(nrow(object$priors$lambda$shape), n_lambda)
  expect_equal(nrow(object$priors$lambda$rate), n_lambda)
  expect_equal(as.numeric(object$priors$lambda$shape[1]), 4)
  expect_equal(as.numeric(object$priors$lambda$rate[1]), 0.05)

  expect_equal(dim(object$priors$a$vinv), rep(n_a, 2))
  expect_equal(nrow(object$priors$a$shape), n_a)
  expect_equal(as.numeric(object$priors$a$rate[1]), 0.06)

  # The error blocks are the gamma ones, unchanged: this model's drift is in the
  # coefficients.
  expect_equal(nrow(object$priors$u$shape), 5)
  expect_equal(nrow(object$priors$v$shape), 2)
})

test_that("add_priors reports a constant specification given to a tvp model", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1, p = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, tvp = TRUE,
                           iterations = 20, burnin = 10)

  # The default arguments carry no state equation, so this is what a caller who
  # forgets gets, and the message has to name what is missing.
  expect_error(add_priors(object), "'shape', 'rate'")
  expect_error(add_priors(object, lambda = tvp_prior(), a = list(vinv = 0.01)),
               "'a' is missing")

  expect_error(add_priors(object, lambda = tvp_prior(rate = 0), a = tvp_prior()),
               "lambda\\$rate")
  expect_error(add_priors(object, lambda = tvp_prior(shape = -1), a = tvp_prior()),
               "lambda\\$shape")

  # A constant-coefficient model is not asked for one.
  constant <- create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 20, burnin = 10)
  expect_silent(add_priors(constant))
})

test_that("add_initial_values produces a path per coefficient block", {

  sim <- sim_dfm(tt = 40, m = 5, n = 2, p = 2)
  object <- create_dfmodel(x = sim$x, p = 2, n = 2, tvp = TRUE,
                           iterations = 20, burnin = 10)
  object <- add_priors(object, lambda = tvp_prior(), a = tvp_prior())
  object <- add_initial_values(object)

  n_lambda <- n_free_lambda(5, 2)
  n_a <- 2 * 2 * 2

  expect_equal(dim(object$initial$lambda), c(n_lambda, sim$tt))
  expect_equal(dim(object$initial$lambda_sigma_inv), rep(n_lambda, 2))
  expect_equal(length(object$initial$lambda_init), n_lambda)

  expect_equal(dim(object$initial$a), c(n_a, sim$tt))
  expect_equal(dim(object$initial$a_sigma_inv), rep(n_a, 2))
  expect_equal(length(object$initial$a_init), n_a)

  # The path starts flat at the state it is given, so every column is that state.
  expect_equal(object$initial$lambda[, 1], as.numeric(object$initial$lambda_init))
  expect_equal(object$initial$lambda[, sim$tt], as.numeric(object$initial$lambda_init))

  # A precision, so positive on the diagonal and zero off it.
  expect_true(all(diag(object$initial$lambda_sigma_inv) > 0))
  expect_equal(sum(object$initial$lambda_sigma_inv[lower.tri(object$initial$lambda_sigma_inv)]), 0)

  # And the two error precisions are still the gamma model's.
  expect_equal(dim(object$initial$uinv), c(5, 5))
  expect_equal(dim(object$initial$vinv), c(2, 2))
})

test_that("the posterior carries a coefficient path and a state variance", {

  prepared <- prepared_dfm_tvp(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)
  object <- add_posterior_coefficients(prepared$object)
  sim <- prepared$sim

  expect_null(object$error)
  expect_equal(names(object$posterior),
               c("lambda", "factors", "a", "u_sigma_inv", "v_sigma_inv"))

  # The two paths, tt times wider than the constant model's, with the periods
  # stacked within a row.
  expect_equal(dim(object$posterior$lambda$coeffs), c(20, sim$m * sim$n * sim$tt))
  expect_equal(dim(object$posterior$a$coeffs), c(20, sim$n^2 * sim$p * sim$tt))

  # And the state variance of each, one number per coefficient.
  expect_equal(dim(object$posterior$lambda$sigma), c(20, n_free_lambda(sim$m, sim$n)))
  expect_equal(dim(object$posterior$a$sigma), c(20, sim$n^2 * sim$p))
  expect_true(all(object$posterior$lambda$sigma > 0))

  # Everything the constant model returns is still the width it was.
  expect_equal(dim(object$posterior$factors$coeffs), c(20, sim$n * sim$tt))
  expect_equal(dim(object$posterior$u_sigma_inv$coeffs), c(20, sim$m))
  expect_equal(dim(object$posterior$v_sigma_inv$coeffs), c(20, sim$n))

  # Both elements of both blocks are chains.
  for (block in c("lambda", "a")) {
    for (element in c("coeffs", "sigma")) {
      expect_s3_class(object$posterior[[block]][[element]], "mcmc")
    }
  }

  # The identifying block is not drawn, so it is exactly itself in every period
  # of every draw and not merely on average.
  path <- matrix(object$posterior$lambda$coeffs[1, ], sim$m * sim$n, sim$tt)
  expect_true(all(path[1, ] == 1))
})

test_that("a transition of order zero leaves the transition block out", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1, p = 1)
  object <- create_dfmodel(x = sim$x, p = 0, n = 1, tvp = TRUE,
                           iterations = 20, burnin = 10)

  # Nothing to give a state equation to, so none is asked for.
  object <- add_priors(object, lambda = tvp_prior())
  expect_null(object$priors$a)

  object <- add_initial_values(object)
  expect_null(object$initial$a)

  object <- add_posterior_coefficients(object)
  expect_null(object$error)
  expect_null(object$posterior$a)

  # The loadings still drift, which is the whole of the model that is left.
  expect_equal(dim(object$posterior$lambda$coeffs), c(20, 4 * 1 * 40))
})

test_that("forecasts and the log likelihood take their own slice of the path", {

  prepared <- prepared_dfm_tvp(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)
  sim <- prepared$sim

  object <- add_posterior_coefficients(prepared$object)
  object <- add_posterior_forecasts(object, n_ahead = 5)
  object <- add_posterior_loglik(object)

  # The forecast holds the coefficients at their last in-sample period, so the
  # binding cuts the path to it; the shape is the constant model's.
  expect_equal(dim(object$posterior$forecast), c(20, 5 * sim$m))
  expect_true(all(is.finite(object$posterior$forecast)))
  expect_s3_class(object$posterior$forecast, "mcmc")

  # The log likelihood scores every period under its own loadings and takes the
  # whole path.
  expect_equal(dim(object$posterior$loglik), c(20, sim$tt))
  expect_true(all(is.finite(object$posterior$loglik)))

  # And neither runs without the factor path, which is part of this posterior
  # rather than derivable from it.
  without <- object
  without$posterior$factors <- NULL
  expect_error(add_posterior_forecasts(without, n_ahead = 5), "factors")
  expect_error(add_posterior_loglik(without), "factors")
})

test_that("it finds loadings that moved over the sample", {

  sim <- sim_dfm_drifting(tt = 400)

  object <- create_dfmodel(x = sim$x, p = 1, n = 1, normalize_x = FALSE,
                           tvp = TRUE, iterations = 300, burnin = 200)
  object <- add_priors(object,
                       lambda = tvp_prior(vinv = 0.01, shape = 3, rate = 0.01),
                       a = tvp_prior(vinv = 0.01, shape = 3, rate = 0.01),
                       u = list(shape = 3, rate = 0.5),
                       v = list(shape = 3, rate = 0.5))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  path <- loading_path(object, sim$m, sim$n, sim$tt) # M x T at one factor
  quarter <- sim$tt %/% 4
  early <- rowMeans(path[, seq_len(quarter), drop = FALSE])
  late <- rowMeans(path[, sim$tt - seq_len(quarter) + 1, drop = FALSE])

  # The level, over the whole sample, against the average of the true ramp.
  expect_equal(rowMeans(path), colMeans(sim$lambda), tolerance = 0.2)

  # The direction of travel, which is what a constant-loading model cannot get
  # right and what this model exists for. Series 2 and 4 fall, series 5 rises,
  # series 1 and 3 do not move.
  expect_lt(late[2], early[2] - 0.3)
  expect_lt(late[4], early[4] - 0.3)
  expect_gt(late[5], early[5] + 0.3)
  expect_equal(late[3], early[3], tolerance = 0.25)

  # The first loading is the identification and does not move at all.
  expect_true(all(path[1, ] == 1))

  # The idiosyncratic variance, which does not drift in this model and should
  # come back at what generated the sample.
  expect_equal(as.numeric(1 / colMeans(object$posterior$u_sigma_inv$coeffs)),
               rep(sim$u_sd^2, sim$m), tolerance = 0.3)
})
