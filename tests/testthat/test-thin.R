# Thinning in the sampler. Which draws are kept is the only thing it may change,
# so a chain run with thin = 3 has to be exactly every third draw of the chain
# run from the same seed and the same starting values with three times the kept
# draws.

# Every block of draws of a posterior, flattened to "block$element" names.
thin_blocks <- function(posterior, prefix = NULL) {
  blocks <- list()
  for (name in names(posterior)) {
    element <- posterior[[name]]
    label <- paste(c(prefix, name), collapse = "$")
    if (is.list(element) && !inherits(element, "mcmc")) {
      blocks <- c(blocks, thin_blocks(element, label))
    } else if (!is.null(element)) {
      blocks[[label]] <- element
    }
  }
  blocks
}

with_sampler_thin <- function(object, iterations, thin) {
  object$model$iterations <- as.integer(iterations)
  object$model$thin <- as.integer(thin)
  object
}

expect_thinned_chain <- function(base, draws = 30, thin = 3) {
  # Built before the seed is set. Left as a promise, `base` would be built on its
  # first use, after set.seed() -- and building it draws starting values, so the
  # unthinned chain would start from a different generator state than the
  # thinned one.
  force(base)
  set.seed(20260913)
  full <- add_posterior_coefficients(base)
  set.seed(20260913)
  thinned <- add_posterior_coefficients(with_sampler_thin(base, draws / thin, thin))

  keep <- seq(thin, draws, thin)
  before <- thin_blocks(full$posterior)
  after <- thin_blocks(thinned$posterior)
  expect_named(after, names(before))
  for (name in names(before)) {
    expect_identical(c(as.matrix(after[[name]])),
                     c(as.matrix(before[[name]])[keep, , drop = FALSE]), label = name)
    expect_equal(attr(after[[name]], "mcpar"), c(thin, draws, thin), label = name)
  }
  thinned
}

test_that("create_dfmodel() and create_favarmodel() take the thinning interval", {
  sim <- sim_dfm()
  expect_null(create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 10, burnin = 5)$model$thin)
  expect_null(create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 10, burnin = 5,
                             thin = 1)$model$thin)
  expect_identical(create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 10, burnin = 5,
                                  thin = 3)$model$thin, 3L)
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1, thin = 0), "single positive integer")

  favar <- make_favar_sample()
  expect_identical(create_favarmodel(x = favar$x, y = favar$yts, p = 1, n = 1,
                                     iterations = 10, burnin = 5, thin = 2)$model$thin, 2L)
  expect_error(create_favarmodel(x = favar$x, y = favar$yts, p = 1, n = 1, thin = 1.5),
               "single positive integer")
})

test_that("a thinned dynamic factor chain is every thin-th draw of the unthinned one", {
  expect_thinned_chain(prepared_dfm(iterations = 30, burnin = 10)$object)
  thinned <- expect_thinned_chain(prepared_dfm_tvp_sv(iterations = 30, burnin = 10)$object)

  # The forecast takes its labels from the model as the coefficients do.
  set.seed(1)
  forecast <- add_posterior_forecasts(thinned, n_ahead = 2)$posterior$forecast
  expect_equal(attr(forecast, "mcpar"), c(3, 30, 3))
})

test_that("a thinned FAVAR chain is every thin-th draw of the unthinned one", {
  favar <- make_favar_sample(tt = 60, n_x = 5)
  model <- create_favarmodel(x = favar$x, y = favar$yts, p = 1, n = 1,
                             normalize_x = FALSE, iterations = 30, burnin = 10)
  model <- add_initial_values(add_priors(model))
  expect_thinned_chain(model)
})
