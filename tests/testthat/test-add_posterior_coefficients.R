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

test_that("the sampler draws from R's RNG, so set.seed reproduces the draws", {

  # This is the test that guards the Armadillo wiring: the vendored core reaches
  # Armadillo through bayests/arma.h, which src/Makevars points at
  # RcppArmadillo, and that is what puts Armadillo's RNG on R's. Losing it
  # compiles, links and runs -- and silently stops honouring set.seed().
  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 2, p = 1)

  set.seed(11)
  first <- add_posterior_coefficients(prep$object)$posterior

  set.seed(11)
  second <- add_posterior_coefficients(prep$object)$posterior

  expect_equal(first, second)

  set.seed(12)
  third <- add_posterior_coefficients(prep$object)$posterior

  expect_false(isTRUE(all.equal(as.numeric(first$factors$coeffs),
                               as.numeric(third$factors$coeffs))))
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
