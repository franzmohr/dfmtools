test_that("add_posterior_coefficients returns one mcmc block per parameter", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 2, p = 2)

  set.seed(1)
  object <- add_posterior_coefficients(prep$object)

  expect_s3_class(object, "dfmodel")
  expect_false(isTRUE(object$error))
  expect_named(object$posterior,
               c("lambda", "factors", "a", "u_sigma_inv", "v_sigma_inv"))

  m <- 4
  n <- 2
  p <- 2
  tt <- 40
  iterations <- 20

  for (i in names(object$posterior)) {
    expect_s3_class(object$posterior[[i]]$coeffs, "mcmc")
    expect_equal(nrow(object$posterior[[i]]$coeffs), iterations)
  }

  # The whole M x N loading matrix, the fixed elements included, so that a draw
  # is reshaped rather than unpacked.
  expect_equal(ncol(object$posterior$lambda$coeffs), m * n)

  # A whole T-period path of every factor per draw.
  expect_equal(ncol(object$posterior$factors$coeffs), tt * n)

  expect_equal(ncol(object$posterior$a$coeffs), n * n * p)
  expect_equal(ncol(object$posterior$u_sigma_inv$coeffs), m)
  expect_equal(ncol(object$posterior$v_sigma_inv$coeffs), n)

  expect_true(all(is.finite(object$posterior$factors$coeffs)))
  expect_true(all(is.finite(object$posterior$lambda$coeffs)))
})

test_that("the identifying block of the loading matrix is fixed in every draw", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 5, n = 3, p = 1)

  set.seed(2)
  object <- add_posterior_coefficients(prep$object)

  m <- 5
  n <- 3

  # Only the product lambda %*% f_t is identified, so the leading N x N block is
  # unit lower triangular and is not drawn.
  for (i in seq_len(nrow(object$posterior$lambda$coeffs))) {
    draw <- matrix(object$posterior$lambda$coeffs[i, ], m, n)
    block <- draw[seq_len(n), , drop = FALSE]
    expect_equal(diag(block), rep(1, n))
    expect_true(all(block[upper.tri(block)] == 0))
  }
})

test_that("the error blocks are precisions and stay positive", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)

  set.seed(3)
  object <- add_posterior_coefficients(prep$object)

  expect_true(all(object$posterior$u_sigma_inv$coeffs > 0))
  expect_true(all(object$posterior$v_sigma_inv$coeffs > 0))
})

test_that("the sampler draws from R's RNG, so a seed reproduces the draws", {

  # This is the test that guards the Armadillo wiring: the vendored core reaches
  # Armadillo through bayests/arma.h, which src/Makevars points at
  # RcppArmadillo, and that is what puts Armadillo's RNG on R's. Losing it
  # compiles, links and runs -- and silently stops honouring the model's seed
  # and set.seed() alike.
  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 2, p = 1)

  first <- add_posterior_coefficients(add_seed(prep$object, 11))$posterior
  second <- add_posterior_coefficients(add_seed(prep$object, 11))$posterior
  expect_equal(first, second)

  third <- add_posterior_coefficients(add_seed(prep$object, 12))$posterior
  expect_false(isTRUE(all.equal(as.numeric(first$factors$coeffs),
                               as.numeric(third$factors$coeffs))))

  # Without a seed of its own the model draws from R's generator as it stands,
  # which is the wiring with nothing set on top of it.
  unseeded <- prep$object
  unseeded$model$seed <- NULL

  set.seed(11)
  fourth <- add_posterior_coefficients(unseeded)$posterior
  set.seed(11)
  fifth <- add_posterior_coefficients(unseeded)$posterior
  expect_equal(fourth, fifth)

  set.seed(12)
  sixth <- add_posterior_coefficients(unseeded)$posterior
  expect_false(isTRUE(all.equal(as.numeric(fourth$factors$coeffs),
                               as.numeric(sixth$factors$coeffs))))
})

test_that("a failed simulation returns the model with error = TRUE", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)

  object <- suppressWarnings(
    add_posterior_coefficients(prep$object,
                               posterior_function = function(x) stop("no draws"))
  )

  expect_s3_class(object, "dfmodel")
  expect_true(object$error)
  expect_null(object$posterior)
})

test_that("a missing algorithm is reported rather than assumed", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)

  object <- prep$object
  object$model$algorithm <- NULL

  expect_true(suppressWarnings(add_posterior_coefficients(object))$error)
})

test_that("a user-supplied posterior_function is used instead of the internal one", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)

  object <- add_posterior_coefficients(prep$object, posterior_function = function(x) {
    x$posterior <- list(marker = TRUE)
    x
  })

  expect_true(object$posterior$marker)
})

test_that("the sampler reports its progress only when asked", {

  # The reporter behind this had been in the package since the C++ core
  # arrived, with nothing to turn it on: every binding built it silent.
  object <- prepared_dfm(iterations = 40, burnin = 20)$object

  quiet <- capture.output(silent <- add_posterior_coefficients(object))
  loud <- capture.output(reported <- add_posterior_coefficients(object, verbose = TRUE))

  expect_length(quiet, 0L)
  expect_true(any(grepl("Progress", loud)))
  expect_true(any(grepl("100%", loud)))

  # Reporting is not sampling: the same model gives the same draws either way.
  expect_identical(as.matrix(silent$posterior$lambda$coeffs),
                   as.matrix(reported$posterior$lambda$coeffs))

  # A model says whether it was watched, and a quiet one carries no field.
  expect_true(reported$model$verbose)
  expect_null(silent$model$verbose)

  expect_error(add_posterior_coefficients(object, verbose = "yes"),
               "must be TRUE or FALSE")
  expect_error(add_posterior_coefficients(object, verbose = c(TRUE, TRUE)),
               "must be TRUE or FALSE")
})

test_that("a factor augmented VAR reports its progress the same way", {

  sim <- make_favar_sample(tt = 60, n_x = 5)
  model <- add_initial_values(add_priors(
    create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1, normalize_x = FALSE,
                      iterations = 40, burnin = 20)))

  # Assigned inside capture.output(), or the returned model is printed and
  # counted as output.
  expect_length(capture.output(quiet <- add_posterior_coefficients(model)), 0L)
  expect_true(any(grepl("Progress",
                        capture.output(loud <- add_posterior_coefficients(model, verbose = TRUE)))))
  expect_identical(as.matrix(quiet$posterior$lambda$coeffs),
                   as.matrix(loud$posterior$lambda$coeffs))
})
