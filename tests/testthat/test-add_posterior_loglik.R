test_that("add_posterior_loglik returns draws by observations", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 2, p = 1)

  set.seed(31)
  object <- add_posterior_loglik(add_posterior_coefficients(prep$object))

  expect_s3_class(object, "dfmodel")

  # The layout waic and loo expect: draws along the rows, observations along the
  # columns.
  expect_equal(dim(object$posterior$loglik), c(20L, 40L))
  expect_true(all(is.finite(object$posterior$loglik)))

  # A density of a 4-variable observation: negative for data this diffuse, and
  # in any case never NA.
  expect_false(anyNA(object$posterior$loglik))
})

test_that("the log-likelihood is the measurement density at the drawn factors", {

  # Two factors, so that the layout of the reshapes is actually tested: at one
  # factor both orderings coincide.
  prep <- prepared_dfm(iterations = 10, burnin = 5, tt = 40, m = 4, n = 2, p = 1)

  set.seed(32)
  object <- add_posterior_loglik(add_posterior_coefficients(prep$object))

  m <- 4
  n <- 2
  tt <- 40
  x <- object$data$x

  # Recomputed in R for the first draw. It is conditional on the factor path, so
  # nothing has to be filtered: lambda, U and f_t are all in the draw.
  #
  # lambda is the M x N matrix column by column; factors is the path period by
  # period, all N factors of a period together.
  draw <- 1
  lambda <- matrix(object$posterior$lambda$coeffs[draw, ], m, n)
  f <- t(matrix(object$posterior$factors$coeffs[draw, ], n, tt))
  u_sigma <- 1 / object$posterior$u_sigma_inv$coeffs[draw, ]

  expected <- vapply(seq_len(tt), function(i) {
    sum(stats::dnorm(x[i, ], mean = lambda %*% f[i, ], sd = sqrt(u_sigma), log = TRUE))
  }, numeric(1))

  expect_equal(as.numeric(object$posterior$loglik[draw, ]), expected,
               tolerance = 1e-8)
})

test_that("add_posterior_loglik rejects a model without draws", {

  prep <- prepared_dfm(iterations = 20, burnin = 10, tt = 40, m = 4, n = 1, p = 1)

  expect_error(add_posterior_loglik(prep$object), "does not contain posterior draws")

  set.seed(33)
  drawn <- add_posterior_coefficients(prep$object)

  without_factors <- drawn
  without_factors$posterior$factors <- NULL
  expect_error(add_posterior_loglik(without_factors),
               "does not contain posterior draws of the factors")
})
