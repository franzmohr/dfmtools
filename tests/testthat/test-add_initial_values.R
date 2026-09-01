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

test_that("starting values are drawn from R's RNG, so set.seed fixes them", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- add_priors(create_dfmodel(x = sim$x, p = 1, n = 1))

  set.seed(7)
  first <- add_initial_values(object)$initial

  set.seed(7)
  second <- add_initial_values(object)$initial

  expect_equal(first, second)
})
