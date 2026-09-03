# The model in which everything moves, end to end through the R interface.
#
# What the numerics do is covered upstream, in BayesTS's
# test/unit_dfm_tvp_stochvol.cpp. What is left for here is the translation, and
# this model adds nothing to the file format that its two parents did not each
# bring -- which is exactly why it is worth a test of its own on this side. The
# two specifications have to *compose*:
#
#   * `create_dfmodel(error = "sv", tvp = TRUE)` has to reach a fourth algorithm
#     rather than one of the three that came before it.
#   * `add_priors()` has to put a state equation on each coefficient block and a
#     volatility block on each error term, in the same call, without either
#     branch treading on the other. All four groups then carry a `shape` and a
#     `rate`, at four different widths, and each pair belongs to the random walk
#     of the block it sits in.
#   * `add_initial_values()` has to produce four paths.
#   * the posterior widens on both sides at once: `lambda` and `a` by tt because
#     the coefficients drift, `u_sigma_inv` and `v_sigma_inv` by tt because the
#     variances do.
#
# One recovery test at the end, for the claim that carrying both drifts is not
# redundant: the sample it runs on has loadings that fall and an idiosyncratic
# volatility that rises, and the model has to put each where it belongs.

test_that("create_dfmodel reaches the fourth algorithm", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1, p = 1)

  object <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv", tvp = TRUE,
                           iterations = 20, burnin = 10)

  expect_s3_class(object, "dfmodel")
  expect_true(object$model$tvp)
  expect_equal(object$model$error, "sv")
  expect_equal(object$model$algorithm, "DfmTvpStochvol")

  # And the other three are where they were: the two arguments are independent.
  spec <- function(error, tvp) {
    create_dfmodel(x = sim$x, p = 1, n = 1, error = error, tvp = tvp,
                   iterations = 20, burnin = 10)$model$algorithm
  }
  expect_equal(spec("gamma", FALSE), "DfmNormalGamma")
  expect_equal(spec("sv", FALSE), "DfmNormalStochvol")
  expect_equal(spec("gamma", TRUE), "DfmTvpGamma")
})

test_that("add_priors builds all four random walks at their own widths", {

  sim <- sim_dfm(tt = 40, m = 5, n = 2, p = 2)
  object <- create_dfmodel(x = sim$x, p = 2, n = 2, error = "sv", tvp = TRUE,
                           iterations = 20, burnin = 10)

  object <- add_priors(object,
                       lambda = tvp_prior(shape = 4, rate = 0.05),
                       a = tvp_prior(shape = 5, rate = 0.06),
                       u = sv_prior(shape = 6, rate = 0.07),
                       v = sv_prior(shape = 7, rate = 0.08))

  n_lambda <- n_free_lambda(5, 2)
  n_a <- 2 * 2 * 2

  # Four groups, four widths, and a shape/rate pair in each that belongs to that
  # group's own random walk. A branch treading on another would show up here as a
  # pair of the wrong length or the wrong value.
  widths <- list(lambda = n_lambda, a = n_a, u = 5, v = 2)
  shapes <- list(lambda = 4, a = 5, u = 6, v = 7)
  for (group in names(widths)) {
    expect_equal(nrow(object$priors[[group]]$shape), widths[[group]], info = group)
    expect_equal(nrow(object$priors[[group]]$rate), widths[[group]], info = group)
    expect_equal(as.numeric(object$priors[[group]]$shape[1]), shapes[[group]], info = group)
  }

  # The coefficient groups keep `vinv`, which is now the precision of the state
  # before the sample; the error groups keep the four fields that are theirs.
  expect_equal(dim(object$priors$lambda$vinv), rep(n_lambda, 2))
  expect_equal(dim(object$priors$a$vinv), rep(n_a, 2))
  for (field in c("mu", "sigma", "offset")) {
    expect_equal(nrow(object$priors$u[[field]]), 5, info = field)
    expect_equal(nrow(object$priors$v[[field]]), 2, info = field)
  }
  expect_null(object$priors$lambda$offset)
  expect_null(object$priors$u$vinv)
})

test_that("add_priors reports whichever half is missing", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1, p = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv", tvp = TRUE,
                           iterations = 20, burnin = 10)

  # Neither default fits this model, and each half has to say so on its own.
  expect_error(add_priors(object), "'shape', 'rate'")
  expect_error(add_priors(object, lambda = tvp_prior(), a = tvp_prior()), "'mu'")
  expect_error(add_priors(object, lambda = tvp_prior(), a = list(vinv = 0.01),
                          u = sv_prior(), v = sv_prior()),
               "'a' is missing")
  expect_error(add_priors(object, lambda = tvp_prior(), a = tvp_prior(),
                          u = sv_prior(), v = list(shape = 5, rate = 4)),
               "'v' is missing")
})

test_that("add_initial_values produces four paths", {

  sim <- sim_dfm(tt = 40, m = 5, n = 2, p = 2)
  prepared <- prepared_dfm_tvp_sv(tt = 40, m = 5, n = 2, p = 2)
  object <- prepared$object

  n_lambda <- n_free_lambda(5, 2)
  n_a <- 2 * 2 * 2

  # The two coefficient paths, one column per period.
  expect_equal(dim(object$initial$lambda), c(n_lambda, sim$tt))
  expect_equal(dim(object$initial$a), c(n_a, sim$tt))
  expect_equal(dim(object$initial$lambda_sigma_inv), rep(n_lambda, 2))
  expect_equal(dim(object$initial$a_sigma_inv), rep(n_a, 2))

  # The two volatility paths, one *row* per period -- the other orientation, and
  # the one the mixture routine works in.
  expect_equal(dim(object$initial$u_h), c(sim$tt, 5))
  expect_equal(dim(object$initial$v_h), c(sim$tt, 2))
  expect_equal(length(object$initial$u_h_init), 5)
  expect_equal(length(object$initial$v_h_init), 2)

  # And no trace of the precisions a constant-variance model would carry.
  expect_null(object$initial$uinv)
  expect_null(object$initial$vinv)
})

test_that("the posterior widens on both sides at once", {

  prepared <- prepared_dfm_tvp_sv(iterations = 20, burnin = 10,
                                  tt = 40, m = 4, n = 1, p = 1)
  object <- add_posterior_coefficients(prepared$object)
  sim <- prepared$sim

  expect_null(object$error)

  # The coefficients are paths because they drift; the precisions are paths
  # because the variances do. This is the only model in which all four are.
  expect_equal(dim(object$posterior$lambda$coeffs), c(20, sim$m * sim$n * sim$tt))
  expect_equal(dim(object$posterior$a$coeffs), c(20, sim$n^2 * sim$p * sim$tt))
  expect_equal(dim(object$posterior$u_sigma_inv$coeffs), c(20, sim$m * sim$tt))
  expect_equal(dim(object$posterior$v_sigma_inv$coeffs), c(20, sim$n * sim$tt))

  # The two state variances stay points.
  expect_equal(dim(object$posterior$lambda$sigma), c(20, n_free_lambda(sim$m, sim$n)))
  expect_equal(dim(object$posterior$a$sigma), c(20, sim$n^2 * sim$p))
  expect_true(all(object$posterior$lambda$sigma > 0))

  expect_equal(dim(object$posterior$factors$coeffs), c(20, sim$n * sim$tt))

  for (block in c("lambda", "a")) {
    for (element in c("coeffs", "sigma")) {
      expect_s3_class(object$posterior[[block]][[element]], "mcmc")
    }
  }

  # The identifying block is not drawn, in any period of any draw.
  path <- matrix(object$posterior$lambda$coeffs[1, ], sim$m * sim$n, sim$tt)
  expect_true(all(path[1, ] == 1))
})

test_that("forecasts and the log likelihood run over the wider posterior", {

  prepared <- prepared_dfm_tvp_sv(iterations = 20, burnin = 10,
                                  tt = 40, m = 4, n = 1, p = 1)
  sim <- prepared$sim

  object <- add_posterior_coefficients(prepared$object)
  object <- add_posterior_forecasts(object, n_ahead = 5)
  object <- add_posterior_loglik(object)

  # Everything is held at its last in-sample period, so the forecast is the same
  # shape a constant model's is.
  expect_equal(dim(object$posterior$forecast), c(20, 5 * sim$m))
  expect_true(all(is.finite(object$posterior$forecast)))
  expect_s3_class(object$posterior$forecast, "mcmc")

  # The log likelihood scores every period under its own loadings and its own
  # precision, and there is one column per period either way.
  expect_equal(dim(object$posterior$loglik), c(20, sim$tt))
  expect_true(all(is.finite(object$posterior$loglik)))

  without <- object
  without$posterior$factors <- NULL
  expect_error(add_posterior_forecasts(without, n_ahead = 5), "factors")
  expect_error(add_posterior_loglik(without), "factors")
})

test_that("it tells a loading that fell from a variance that rose", {

  # The loadings ramp down and the idiosyncratic scale ramps up, which pull the
  # same way on how large an observed series is and opposite ways on how much of
  # it the factor explains. A model carrying only one of the two drifts would
  # have to explain the other with what it has.
  sim <- sim_dfm_drifting(tt = 400, u_sd = 0.3, u_sd_to = 0.9)

  object <- create_dfmodel(x = sim$x, p = 1, n = 1, normalize_x = FALSE,
                           error = "sv", tvp = TRUE,
                           iterations = 300, burnin = 200)
  object <- add_priors(object,
                       lambda = tvp_prior(vinv = 0.01, shape = 3, rate = 0.01),
                       a = tvp_prior(vinv = 0.01, shape = 3, rate = 0.01),
                       u = sv_prior(v_i = 0.1, shape = 3, rate = 0.01),
                       v = sv_prior(v_i = 0.1, shape = 3, rate = 0.01))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  quarter <- sim$tt %/% 4
  early <- seq_len(quarter)
  late <- sim$tt - seq_len(quarter) + 1

  # The loadings, which fell.
  path <- loading_path(object, sim$m, sim$n, sim$tt)
  lambda_early <- rowMeans(path[, early, drop = FALSE])
  lambda_late <- rowMeans(path[, late, drop = FALSE])

  expect_lt(lambda_late[2], lambda_early[2] - 0.3)
  expect_lt(lambda_late[4], lambda_early[4] - 0.3)
  expect_gt(lambda_late[5], lambda_early[5] + 0.3)
  expect_true(all(path[1, ] == 1))

  # The idiosyncratic variance, which rose -- at the same time, on the same data.
  variance <- idiosyncratic_variance(object, sim$m, sim$tt)
  u_early <- mean(variance[, early])
  u_late <- mean(variance[, late])

  expect_gt(u_late, u_early)
  expect_equal(u_early, sim$u_sd^2, tolerance = 0.6)
  expect_equal(u_late, sim$u_sd_to^2, tolerance = 0.5)
})
