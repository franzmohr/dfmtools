test_that("add_priors sizes every prior from the model dimensions", {

  sim <- sim_dfm(tt = 40, m = 4, n = 2)
  object <- create_dfmodel(x = sim$x, p = 2, n = 2)
  object <- add_priors(object)

  m <- 4
  n <- 2
  p <- 2

  # The loading prior covers the free elements only, ordered row by row, which is
  # the order the equation-by-equation draw consumes them in.
  n_lambda <- n_free_lambda(m, n)
  expect_equal(n_lambda, 5)
  expect_equal(dim(object$priors$lambda$vinv), c(n_lambda, n_lambda))
  expect_equal(diag(object$priors$lambda$vinv), rep(0.01, n_lambda))

  expect_equal(dim(object$priors$a$mu), c(n * n * p, 1L))
  expect_equal(dim(object$priors$a$vinv), c(n * n * p, n * n * p))
  expect_true(all(object$priors$a$mu == 0))

  expect_equal(dim(object$priors$u$shape), c(m, 1L))
  expect_equal(dim(object$priors$u$rate), c(m, 1L))
  expect_equal(dim(object$priors$v$shape), c(n, 1L))
  expect_equal(dim(object$priors$v$rate), c(n, 1L))
})

test_that("add_priors passes the specified values through", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1)
  object <- add_priors(object,
                       lambda = list(vinv = 0.5),
                       u = list(shape = 3, rate = 2),
                       a = list(vinv = 0.25),
                       v = list(shape = 7, rate = 6))

  expect_equal(diag(object$priors$lambda$vinv), rep(0.5, n_free_lambda(4, 1)))
  expect_equal(diag(object$priors$a$vinv), 0.25)
  expect_true(all(object$priors$u$shape == 3))
  expect_true(all(object$priors$u$rate == 2))
  expect_true(all(object$priors$v$shape == 7))
  expect_true(all(object$priors$v$rate == 6))
})

test_that("a transition of order zero gets no coefficient prior", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 0, n = 1)
  object <- add_priors(object)

  expect_null(object$priors$a)
  expect_false(is.null(object$priors$lambda))
})

test_that("add_priors dispatches over a modellist", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- add_priors(create_dfmodel(x = sim$x, p = 1:2, n = 1))

  expect_s3_class(object, "modellist")
  expect_length(object, 2)
  expect_true(all(vapply(object, function(z) !is.null(z$priors), logical(1))))

  # The transition prior tracks each model's own p.
  expect_equal(vapply(object, function(z) nrow(z$priors$a$vinv), numeric(1)),
               c(1, 2))
})

test_that("add_priors rejects invalid input", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1)

  expect_error(add_priors(object, lambda = list()), "'lambda\\$vinv' is missing")
  expect_error(add_priors(object, lambda = list(vinv = -1)), "must be at least 0")
  expect_error(add_priors(object, a = list()), "'a\\$vinv' is missing")
  expect_error(add_priors(object, a = list(vinv = -1)), "must be at least 0")

  expect_error(add_priors(object, u = list(shape = 5)), "at least of length 2")
  expect_error(add_priors(object, u = list(shape = 5, scale = 4)), "u\\$rate is missing")
  expect_error(add_priors(object, u = list(shape = -1, rate = 4)), "'u\\$shape' must be at least 0")
  expect_error(add_priors(object, u = list(shape = 5, rate = 0)), "'u\\$rate' must be larger than 0")

  expect_error(add_priors(object, v = list(shape = 5)), "at least of length 2")
  expect_error(add_priors(object, v = list(shape = 5, scale = 4)), "v\\$rate is missing")
  expect_error(add_priors(object, v = list(shape = -1, rate = 4)), "'v\\$shape' must be at least 0")
  expect_error(add_priors(object, v = list(shape = 5, rate = 0)), "'v\\$rate' must be larger than 0")
})
