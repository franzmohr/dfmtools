test_that("add_initial_values produces one starting value per drawn block", {

  sim <- sim_dfm(tt = 40, m = 4, n = 2)
  object <- add_initial_values(add_priors(create_dfmodel(x = sim$x, p = 2, n = 2)))

  m <- 4
  n <- 2
  p <- 2

  expect_named(object$initial, c("lambda", "uinv", "vinv", "a"))

  # The free loadings only, in the row-by-row order the draw consumes them in.
  expect_equal(dim(object$initial$lambda), c(n_free_lambda(m, n), 1L))
  expect_equal(dim(object$initial$a), c(n * n * p, 1L))

  # Precisions, and diagonal: only the diagonal of either covariance matrix is
  # estimated under the gamma specification.
  expect_equal(dim(object$initial$uinv), c(m, m))
  expect_equal(dim(object$initial$vinv), c(n, n))
  expect_true(all(object$initial$uinv[upper.tri(object$initial$uinv)] == 0))
  expect_true(all(object$initial$vinv[upper.tri(object$initial$vinv)] == 0))
  expect_true(all(diag(object$initial$uinv) > 0))
  expect_true(all(diag(object$initial$vinv) > 0))
})

test_that("a transition of order zero gets no starting value for a", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- add_initial_values(add_priors(create_dfmodel(x = sim$x, p = 0, n = 1)))

  expect_null(object$initial$a)
  expect_false(is.null(object$initial$lambda))
})

test_that("starting values are drawn at the scale of the prior, not its inverse", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  # Every prior here is a precision, so a draw from it has standard deviation
  # 1/sqrt(vinv) -- which is what backsolve(chol(vinv), z) gives. The other
  # arrangement of the same two objects, chol(vinv) %*% z, gives sqrt(vinv), and
  # this file held both at once until they were reconciled.
  #
  # The two agree at vinv = 1 and nowhere else, and they move in opposite
  # directions as the prior tightens, so the check runs at a diffuse precision
  # and at a tight one. A tight prior is the case that matters: there the wrong
  # form throws the starting values far outside the prior rather than merely
  # bunching them inside it.
  scale_of <- function(vinv, reps = 1000) {
    object <- add_priors(create_dfmodel(x = sim$x, p = 1, n = 1),
                         lambda = list(vinv = vinv), a = list(vinv = vinv))
    draws <- replicate(reps, {
      initial <- add_initial_values(object)$initial
      c(initial$lambda[1], initial$a[1])
    })
    apply(draws, 1, stats::sd)
  }

  set.seed(20)
  expect_equal(scale_of(0.01), c(10, 10), tolerance = 0.1)
  expect_equal(scale_of(100), c(0.1, 0.1), tolerance = 0.1)
})

test_that("starting values are drawn from R's RNG, so set.seed fixes them", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- add_priors(create_dfmodel(x = sim$x, p = 1, n = 1))

  set.seed(7)
  first <- add_initial_values(object)$initial

  set.seed(7)
  second <- add_initial_values(object)$initial

  expect_equal(first, second)
})

# One observed series and one factor is the only specification whose loading
# matrix is entirely fixed: lambda is 1 x 1, the identifying restriction puts a
# one in it, and (2M - N - 1)N/2 is zero. What is left is x_t = f_t + u_t with an
# AR(p) factor, which is an unobserved-components model and which the sampler
# estimates like any other -- so the empty block is carried through rather than
# rejected. Before this was handled, chol() met the 0 x 0 prior precision and
# reported a zero-dimensional matrix instead of anything about the model.
single_series <- function(tt = 60) {
  set.seed(11)
  f <- numeric(tt + 20)
  for (i in 2:(tt + 20)) f[i] <- 0.7 * f[i - 1] + stats::rnorm(1)
  f <- f[-(1:20)]
  stats::ts(matrix(f + stats::rnorm(tt, sd = 0.5), ncol = 1,
                   dimnames = list(NULL, "a")),
            start = c(1980, 1), frequency = 4)
}

test_that("one series and one factor leaves an empty loading block, not an error", {

  x <- single_series()

  object <- add_priors(create_dfmodel(x = x, p = 1, n = 1,
                                      iterations = 20, burnin = 10))
  expect_equal(n_free_lambda(1, 1), 0)
  expect_equal(dim(object$priors$lambda$vinv), c(0L, 0L))

  object <- add_initial_values(object)

  expect_equal(dim(object$initial$lambda), c(0L, 1L))
  expect_equal(dim(object$initial$a), c(1L, 1L))
  expect_equal(dim(object$initial$uinv), c(1L, 1L))
  expect_equal(dim(object$initial$vinv), c(1L, 1L))

  # And the sampler runs it, which is why the specification is supported rather
  # than turned away. The whole of lambda comes back, as it does for every other
  # model here, and every draw of it is the fixed one.
  object <- add_posterior_coefficients(object)
  expect_equal(dim(object$posterior$lambda$coeffs), c(20L, 1L))
  expect_true(all(object$posterior$lambda$coeffs == 1))
  expect_equal(dim(object$posterior$a$coeffs), c(20L, 1L))
})

test_that("the empty loading block is empty for time varying coefficients too", {

  x <- single_series()
  tt <- nrow(x)

  object <- create_dfmodel(x = x, p = 1, n = 1, tvp = TRUE,
                           iterations = 20, burnin = 10)
  object <- add_priors(object, lambda = tvp_prior(), a = tvp_prior())

  # The state equation of the loadings has nothing to describe either.
  expect_length(object$priors$lambda$shape, 0)
  expect_length(object$priors$lambda$rate, 0)

  object <- add_initial_values(object)

  expect_equal(dim(object$initial$lambda), c(0L, tt))
  expect_equal(dim(object$initial$lambda_init), c(0L, 1L))
  expect_equal(dim(object$initial$lambda_sigma_inv), c(0L, 0L))

  # The transition still drifts, and is the only block that does.
  expect_equal(dim(object$initial$a), c(1L, tt))

  object <- add_posterior_coefficients(object)
  expect_equal(dim(object$posterior$lambda$coeffs), c(20L, tt))
  expect_true(all(object$posterior$lambda$coeffs == 1))
})
