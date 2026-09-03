# A modellist is what create_dfmodel returns for a vector of p or n, i.e. the
# specification search the package is meant to support. The individual tests
# check that add_priors dispatches over such a list; what is untested is the
# rest of the chain, where the failure mode is a list whose elements have
# silently been given each other's dimensions.

test_that("the whole pipeline dispatches over a modellist", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  object <- create_dfmodel(x = sim$x, p = 1:2, n = 1:2,
                           iterations = 20, burnin = 10)
  object <- add_priors(object)
  object <- add_initial_values(object)

  set.seed(61)
  object <- add_posterior_coefficients(object)

  expect_s3_class(object, "modellist")
  expect_length(object, 4)
  expect_false(any(vapply(object, function(z) isTRUE(z$error), logical(1))))

  # Every element keeps its own specification, and its blocks are sized from
  # that specification rather than from the first model in the list.
  for (i in seq_along(object)) {

    model <- object[[i]]
    p <- model$model$p
    n <- model$model$n
    m <- 4
    tt <- 40

    expect_equal(ncol(model$posterior$lambda$coeffs), m * n, info = i)
    expect_equal(ncol(model$posterior$factors$coeffs), tt * n, info = i)
    expect_equal(ncol(model$posterior$a$coeffs), n * n * p, info = i)
    expect_equal(ncol(model$posterior$u_sigma_inv$coeffs), m, info = i)
    expect_equal(ncol(model$posterior$v_sigma_inv$coeffs), n, info = i)
  }

  # And the four elements really are the four combinations, not one repeated.
  specs <- t(vapply(object, function(z) c(z$model$p, z$model$n), numeric(2)))
  expect_equal(nrow(unique(specs)), 4L)
})

test_that("forecasts and the log-likelihood dispatch over a modellist", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  object <- create_dfmodel(x = sim$x, p = 1:2, n = 1,
                           iterations = 20, burnin = 10)
  object <- add_initial_values(add_priors(object))

  set.seed(62)
  object <- add_posterior_coefficients(object)
  object <- add_posterior_forecasts(object, n_ahead = 3)
  object <- add_posterior_loglik(object)

  expect_s3_class(object, "modellist")

  for (i in seq_along(object)) {
    expect_equal(dim(object[[i]]$posterior$forecast), c(20L, 3L * 4L), info = i)
    expect_equal(dim(object[[i]]$posterior$loglik), c(20L, 40L), info = i)
    expect_true(all(is.finite(object[[i]]$posterior$forecast)), info = i)
    expect_true(all(is.finite(object[[i]]$posterior$loglik)), info = i)
  }

  # The two models are fitted to the same data, so their log-likelihoods are
  # comparable -- which is the point of estimating a list of them -- but they
  # are not the same numbers.
  expect_false(isTRUE(all.equal(as.numeric(object[[1]]$posterior$loglik),
                                as.numeric(object[[2]]$posterior$loglik))))
})

test_that("a modellist mixing error specifications keeps them apart", {

  # create_dfmodel produces one error specification per call, so a list mixing
  # the two is assembled by the caller. Each element then has to reach its own
  # sampler, which is dispatched per element from model$algorithm.
  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  gamma_model <- add_initial_values(add_priors(
    create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 20, burnin = 10)))

  sv_model <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv",
                             iterations = 20, burnin = 10)
  sv_model <- add_initial_values(add_priors(sv_model, u = sv_prior(), v = sv_prior()))

  object <- list(gamma_model, sv_model)
  class(object) <- append("modellist", class(object))

  set.seed(63)
  object <- add_posterior_coefficients(object)

  expect_false(any(vapply(object, function(z) isTRUE(z$error), logical(1))))

  # One variance per series against a whole path per series, which is the
  # difference between the two samplers.
  expect_equal(ncol(object[[1]]$posterior$u_sigma_inv$coeffs), 4L)
  expect_equal(ncol(object[[2]]$posterior$u_sigma_inv$coeffs), 4L * 40L)
})
