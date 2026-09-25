# The non-centred prior `omega_v` on the random walks of a time varying model,
# from what add_priors() stores through what the sampler returns to the Bayes
# factors time_variation_test() reads off it.

test_that("add_priors stores omega_v in place of the state equation's gamma", {

  sim <- sim_dfm(tt = 40, m = 4, n = 2)
  object <- create_dfmodel(x = sim$x, p = 2, n = 2, tvp = TRUE)
  object <- add_priors(object,
                       lambda = list(vinv = 0.01, omega_v = 0.1),
                       a = list(vinv = 0.01, omega_v = 0.2))

  n_lambda <- n_free_lambda(4, 2)
  expect_equal(dim(object$priors$lambda$omega_v), c(n_lambda, 1L))
  expect_true(all(object$priors$lambda$omega_v == 0.1))
  expect_equal(dim(object$priors$a$omega_v), c(2L * 2L * 2L, 1L))
  expect_true(all(object$priors$a$omega_v == 0.2))

  # The two parameterisations are alternatives, so the gamma is not there as well.
  expect_null(object$priors$lambda$shape)
  expect_null(object$priors$lambda$rate)
  expect_null(object$priors$a$shape)
  expect_null(object$priors$a$rate)
})

test_that("one block may be non-centred while the other keeps the gamma prior", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, tvp = TRUE)
  object <- add_priors(object,
                       lambda = list(vinv = 0.01, omega_v = 0.1),
                       a = tvp_prior())

  expect_true(all(object$priors$lambda$omega_v == 0.1))
  expect_null(object$priors$a$omega_v)
  expect_true(all(object$priors$a$shape == 3))
})

test_that("a stochastic volatility group takes omega_v and keeps its other elements", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv", tvp = TRUE)
  sv <- sv_prior()
  sv$shape <- NULL
  sv$rate <- NULL
  sv$omega_v <- 0.05
  object <- add_priors(object, lambda = tvp_prior(), a = tvp_prior(), u = sv, v = sv)

  expect_true(all(object$priors$u$omega_v == 0.05))
  expect_null(object$priors$u$shape)
  expect_null(object$priors$u$rate)

  # What the log-volatility block reads either way: the prior on its initial
  # state, the starting value of the innovation variance and the offset.
  expect_equal(dim(object$priors$u$mu), c(4L, 1L))
  expect_equal(diag(object$priors$u$v_inv), rep(1, 4))
  expect_true(all(object$priors$u$sigma == 0.05))
  expect_true(all(object$priors$u$offset == 1e-4))
})

test_that("add_priors refuses omega_v beside the gamma pair and where nothing drifts", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  tvp <- create_dfmodel(x = sim$x, p = 1, n = 1, tvp = TRUE)
  expect_error(
    add_priors(tvp, lambda = list(vinv = 0.01, omega_v = 0.1, shape = 3, rate = 0.01),
               a = tvp_prior()),
    "both 'omega_v' and")
  expect_error(
    add_priors(tvp, lambda = tvp_prior(), a = list(vinv = 0.01, omega_v = 0.1, rate = 0.01)),
    "both 'omega_v' and")

  # A constant-coefficient model has no random walk in either block.
  constant <- create_dfmodel(x = sim$x, p = 1, n = 1)
  expect_error(add_priors(constant, lambda = list(vinv = 0.01, omega_v = 0.1)),
               "only available for a model with tvp = TRUE")
  expect_error(add_priors(constant, a = list(vinv = 0.01, omega_v = 0.1)),
               "only available for a model with tvp = TRUE")

  # Nor does a constant error precision, whatever the coefficients do.
  expect_error(add_priors(tvp, lambda = tvp_prior(), a = tvp_prior(),
                          u = list(shape = 5, rate = 4, omega_v = 0.1)),
               "only available for a model with error = \"sv\" and tvp = TRUE")

  # The log-volatilities of a constant-coefficient model drift, but only
  # DfmTvpStochvol draws one non-centred, so the group is refused here rather
  # than by the core complaining about the shape it is then missing.
  sv <- sv_prior()
  sv$shape <- NULL
  sv$rate <- NULL
  sv$omega_v <- 0.05
  constant_sv <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv")
  expect_error(add_priors(constant_sv, u = sv, v = sv_prior()),
               "only available for a model with error = \"sv\" and tvp = TRUE")
  expect_error(add_priors(constant_sv, u = sv_prior(), v = sv),
               "only available for a model with error = \"sv\" and tvp = TRUE")
})

test_that("omega_v is checked for what it is and where it lies", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, tvp = TRUE)

  expect_error(add_priors(object, lambda = list(vinv = 0.01, omega_v = 0), a = tvp_prior()),
               "must be larger than 0")
  expect_error(add_priors(object, lambda = list(vinv = 0.01, omega_v = -1), a = tvp_prior()),
               "must be larger than 0")
  expect_error(add_priors(object, lambda = list(vinv = 0.01, omega_v = "a"), a = tvp_prior()),
               "must be a single finite number")
})

test_that("add_initial_values draws the state variance from the non-centred prior", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, tvp = TRUE)
  object <- add_priors(object,
                       lambda = list(vinv = 0.01, omega_v = 0.1),
                       a = list(vinv = 0.01, omega_v = 0.1))
  object <- add_initial_values(object)

  n_lambda <- n_free_lambda(4, 1)
  expect_equal(dim(object$initial$lambda_sigma_inv), c(n_lambda, n_lambda))
  expect_true(all(is.finite(diag(object$initial$lambda_sigma_inv))))
  expect_true(all(diag(object$initial$lambda_sigma_inv) > 0))
  expect_equal(dim(object$initial$a_sigma_inv), c(1L, 1L))
  expect_true(all(is.finite(diag(object$initial$a_sigma_inv))))
})

test_that("the sampler returns omega beside sigma, and sigma is its square", {

  sim <- sim_dfm(tt = 40, m = 4, n = 2)
  object <- create_dfmodel(x = sim$x, p = 1, n = 2, tvp = TRUE,
                           iterations = 20, burnin = 10)
  object <- add_priors(object,
                       lambda = list(vinv = 0.01, omega_v = 0.1),
                       a = list(vinv = 0.01, omega_v = 0.1))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  n_lambda <- n_free_lambda(4, 2)
  for (block in c("lambda", "a")) {
    draws <- object$posterior[[block]]
    expect_false(is.null(draws$omega))
    expect_false(is.null(draws$omega_log_zero))
    expect_false(is.null(draws$omega_log_zero_joint))
    expect_equal(ncol(draws$omega_log_zero_joint), 1L)

    # sigma = omega^2 is the parameterisation itself, so the two must agree
    # element by element -- which they can only do if the rows of both went
    # through the same permutation of the free loadings on the way out.
    expect_equal(as.matrix(draws$omega)^2, as.matrix(draws$sigma),
                 ignore_attr = TRUE)
  }
  expect_equal(ncol(object$posterior$lambda$omega), n_lambda)
  expect_equal(ncol(object$posterior$a$omega), 2L * 2L)
})

test_that("a centred block returns no omega", {

  object <- prepared_dfm_tvp(iterations = 20, burnin = 10)$object
  object <- add_posterior_coefficients(object)

  expect_null(object$posterior$lambda$omega)
  expect_null(object$posterior$a$omega_log_zero)
})

test_that("all four random walks of a stochastic volatility model can be non-centred", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  sv <- sv_prior()
  sv$shape <- NULL
  sv$rate <- NULL
  sv$omega_v <- 0.05
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv", tvp = TRUE,
                           iterations = 20, burnin = 10)
  object <- add_priors(object,
                       lambda = list(vinv = 0.01, omega_v = 0.1),
                       a = list(vinv = 0.01, omega_v = 0.1),
                       u = sv, v = sv)
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  for (block in c("lambda", "a", "u_sigma_inv", "v_sigma_inv")) {
    expect_false(is.null(object$posterior[[block]]$omega),
                 info = paste("omega of", block))
  }
  expect_equal(ncol(object$posterior$u_sigma_inv$omega), 4L)
  expect_equal(ncol(object$posterior$v_sigma_inv$omega), 1L)
})

test_that("time_variation_test reports one row per state and a joint row per block", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 2, n = 1, tvp = TRUE,
                           iterations = 40, burnin = 20)
  object <- add_priors(object,
                       lambda = list(vinv = 0.01, omega_v = 0.1),
                       a = list(vinv = 0.01, omega_v = 0.1))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  result <- time_variation_test(object)
  expect_s3_class(result, "bvartimevar")
  expect_equal(colnames(result), c("block", "equation", "term", "log_bf", "nse"))

  n_lambda <- n_free_lambda(4, 1)
  expect_equal(sum(result$block == "loadings"), n_lambda + 1L)
  expect_equal(sum(result$block == "transition"), 2L + 1L)
  expect_equal(sum(result$term == "(joint)"), 2L)
  expect_true(all(is.finite(result$log_bf)))

  # The labels of the loadings are the free elements in R's ordering, which is
  # the column-major order of lower.tri() over the M x N matrix.
  loadings <- as.data.frame(result)[result$block == "loadings" & result$term != "(joint)", ]
  expect_equal(loadings$equation, colnames(sim$x)[2:4])
  expect_equal(unique(loadings$term), "factor1")

  # vec([A_1 A_2]), so the lag is what changes between the two rows here.
  transition <- as.data.frame(result)[result$block == "transition" & result$term != "(joint)", ]
  expect_equal(transition$term, c("factor1.l1", "factor1.l2"))

  without_joint <- time_variation_test(object, joint = FALSE)
  expect_equal(sum(without_joint$term == "(joint)"), 0L)
  expect_equal(nrow(without_joint), nrow(result) - 2L)
})

test_that("time_variation_test refuses a model it cannot test", {

  object <- prepared_dfm_tvp(iterations = 20, burnin = 10)$object
  expect_error(time_variation_test(object), "no posterior draws")

  drawn <- add_posterior_coefficients(object)
  expect_error(time_variation_test(drawn), "No block of the model was estimated")

  non_centred <- create_dfmodel(x = sim_dfm(tt = 40, m = 4, n = 1)$x, p = 1, n = 1,
                                tvp = TRUE, iterations = 20, burnin = 10)
  non_centred <- add_priors(non_centred,
                            lambda = list(vinv = 0.01, omega_v = 0.1),
                            a = list(vinv = 0.01, omega_v = 0.1))
  non_centred <- add_posterior_coefficients(add_initial_values(non_centred))
  expect_error(time_variation_test(non_centred, joint = NA), "must be TRUE or FALSE")
  expect_error(time_variation_test(non_centred, batches = 1), "at least 2")
})

test_that("a block that moves is told from one that does not", {

  # The loadings of this sample really drift and its transition is constant, so
  # the two blocks should not come out the same way. A short chain is enough for
  # the sign, which is all this asserts.
  sim <- sim_dfm_drifting(tt = 200, seed = 11)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, tvp = TRUE,
                           iterations = 200, burnin = 100)
  object <- add_priors(object,
                       lambda = list(vinv = 0.01, omega_v = 0.1),
                       a = list(vinv = 0.01, omega_v = 0.1))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  result <- as.data.frame(time_variation_test(object, joint = FALSE))
  loadings <- result[result$block == "loadings", ]
  transition <- result[result$block == "transition", ]

  expect_true(max(loadings$log_bf) > 3)
  expect_true(max(loadings$log_bf) > max(transition$log_bf))
})
