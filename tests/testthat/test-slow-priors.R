# The slow-moving restriction of add_priors.favarmodel: the loadings a series
# carries on the observed block, held at zero rather than estimated.

# Where series j's loading on observed variable q sits in the free-loading
# vector, which runs row by row over series n + 1 .. k and left to right within
# a row. Written out independently of the function under test.
slow_position <- function(j, q, n, n_state) (j - n - 1) * n_state + n + q

test_that("the restriction lands on the loadings it names and no others", {

  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1,
                             iterations = 10, burnin = 5)

  free <- add_priors(model)
  rest <- add_priors(model, slow = c(3, 5))

  n_state <- 2
  hit <- vapply(c(3, 5), slow_position, numeric(1), q = 1, n = 1,
                n_state = n_state)

  expect_equal(diag(rest$priors$lambda$vinv)[hit], rep(1e12, 2))
  # Everything else is left at the loading prior it would have had.
  expect_equal(diag(rest$priors$lambda$vinv)[-hit],
               diag(free$priors$lambda$vinv)[-hit])
  # The prior mean is zero throughout, so the restriction is to zero.
  expect_true(all(rest$priors$lambda$mu == 0))
})

test_that("series are named as well as indexed", {

  sim <- make_favar_sample()
  colnames(sim$x) <- paste0("v", seq_len(ncol(sim$x)))
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1,
                             iterations = 10, burnin = 5)

  by_name <- add_priors(model, slow = c("v3", "v5"))
  by_index <- add_priors(model, slow = c(3, 5))

  expect_equal(by_name$priors$lambda$vinv, by_index$priors$lambda$vinv)
})

test_that("the identifying series may be listed and change nothing", {

  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 2,
                             iterations = 10, burnin = 5)

  # Series 1 and 2 identify the factors, so they carry no free loading at all.
  with_them <- add_priors(model, slow = c(1, 2, 4))
  without <- add_priors(model, slow = 4)

  expect_equal(with_them$priors$lambda$vinv, without$priors$lambda$vinv)
  # And listing only those is a restriction on nothing.
  expect_equal(add_priors(model, slow = c(1, 2))$priors$lambda$vinv,
               add_priors(model)$priors$lambda$vinv)
})

test_that("every observed column of a slow series is restricted", {

  sim <- make_favar_sample()
  y2 <- stats::ts(cbind(sim$y, rev(sim$y)), start = c(1990, 1), frequency = 4)
  model <- create_favarmodel(x = sim$x, y = y2, p = 1, n = 1,
                             iterations = 10, burnin = 5)

  rest <- add_priors(model, slow = 4)

  n_state <- 3
  hit <- vapply(1:2, slow_position, numeric(1), j = 4, n = 1, n_state = n_state)
  expect_equal(diag(rest$priors$lambda$vinv)[hit], rep(1e12, 2))
  expect_equal(sum(diag(rest$priors$lambda$vinv) == 1e12), 2)
})

test_that("the list form sets the precision the restriction is held at", {

  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1,
                             iterations = 10, burnin = 5)

  rest <- add_priors(model, slow = list(series = 4, vinv = 500))

  expect_equal(diag(rest$priors$lambda$vinv)[slow_position(4, 1, 1, 2)], 500)
  # An empty list is no restriction.
  expect_equal(add_priors(model, slow = list(series = NULL))$priors$lambda$vinv,
               add_priors(model)$priors$lambda$vinv)
})

test_that("the restriction is checked", {

  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1,
                             iterations = 10, burnin = 5)

  expect_error(add_priors(model, slow = "nowhere"),
               "not available as 'slow'")
  expect_error(add_priors(model, slow = 99),
               "between 1 and 8")
  expect_error(add_priors(model, slow = list(series = 4, vinv = -1)),
               "at least 0")
})

test_that("a restricted loading comes back at zero and its neighbours do not", {

  # The point of the restriction is what the sampler then does with it, so it is
  # run rather than only inspected.
  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1,
                             normalize_x = FALSE,
                             iterations = 300, burnin = 100)
  model <- add_priors(model, slow = 2:4)
  model <- add_initial_values(model)
  model <- add_posterior_coefficients(model)

  lambda <- matrix(colMeans(model$posterior$lambda$coeffs), 8, 2)

  # The identifying row is untouched, the restricted rows are at zero on the
  # observed column, and the free ones are not.
  expect_equal(lambda[1, ], c(1, 0))
  expect_true(all(abs(lambda[2:4, 2]) < 1e-6))
  expect_true(all(abs(lambda[5:8, 2]) > 1e-3))
  # The factor loadings of the restricted series stay free.
  expect_true(all(abs(lambda[2:4, 1]) > 1e-3))
})
