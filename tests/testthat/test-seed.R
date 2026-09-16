# The seed of the posterior simulation, as bvartools has it: add_initial_values()
# stores one, add_seed() replaces it, and add_posterior_coefficients() draws with
# it and leaves R's generator as it found it.

prepared_favar <- function(iterations = 20, burnin = 10) {
  favar <- make_favar_sample(tt = 60, n_x = 5)
  model <- create_favarmodel(x = favar$x, y = favar$yts, p = 1, n = 1,
                             normalize_x = FALSE,
                             iterations = iterations, burnin = burnin)
  add_initial_values(add_priors(model))
}

# One model for each of the five samplers.
one_model_per_sampler <- function() {
  list(DfmNormalGamma = prepared_dfm()$object,
       DfmNormalStochvol = prepared_dfm_sv()$object,
       DfmTvpGamma = prepared_dfm_tvp()$object,
       DfmTvpStochvol = prepared_dfm_tvp_sv()$object,
       FavarNormalWishart = prepared_favar())
}

# Models with priors and no starting values, built without touching the
# generator afterwards: sim_dfm() and make_favar_sample() seed it themselves, so
# the set.seed() calls in a test are what decides everything after them.
unseeded_dfm <- function() {
  add_priors(create_dfmodel(x = sim_dfm()$x, p = 1, n = 1,
                            iterations = 20, burnin = 10))
}

unseeded_favar <- function() {
  sample <- make_favar_sample(tt = 60, n_x = 5)
  add_priors(create_favarmodel(x = sample$x, y = sample$yts, p = 1, n = 1,
                               iterations = 20, burnin = 10))
}

test_that("add_initial_values() stores a seed that set.seed() reproduces", {
  for (model in list(unseeded_dfm(), unseeded_favar())) {
    set.seed(1)
    first <- add_initial_values(model)
    set.seed(1)
    again <- add_initial_values(model)
    set.seed(2)
    other <- add_initial_values(model)

    expect_type(first$model$seed, "integer")
    expect_gte(first$model$seed, 0L)
    expect_identical(first$model$seed, again$model$seed)
    expect_false(identical(first$model$seed, other$model$seed))
  }
})

test_that("add_initial_values() keeps a seed the model already has", {
  expect_identical(add_initial_values(add_seed(unseeded_dfm(), 42))$model$seed, 42L)
  expect_identical(add_initial_values(add_seed(unseeded_favar(), 43))$model$seed, 43L)
})

test_that("add_seed() stores a seed as an integer and refuses anything else", {
  model <- prepared_dfm()$object
  expect_identical(add_seed(model, 20260916)$model$seed, 20260916L)
  expect_identical(add_seed(prepared_favar(), 0)$model$seed, 0L)

  for (bad in list(-1, 1.5, NA, "1", c(1, 2), Inf, .Machine$integer.max + 1)) {
    expect_error(add_seed(model, bad), "whole number")
  }
})

test_that("the seed decides the draws of every sampler and leaves R's generator as it was", {
  models <- one_model_per_sampler()

  for (algorithm in names(models)) {
    model <- add_seed(models[[algorithm]], 7)
    expect_identical(model$model$algorithm, algorithm)

    set.seed(1)
    state <- get(".Random.seed", envir = globalenv())
    first <- add_posterior_coefficients(model)
    expect_identical(get(".Random.seed", envir = globalenv()), state)
    expect_null(first$error)

    set.seed(2)
    second <- add_posterior_coefficients(model)
    expect_identical(second$posterior, first$posterior)

    third <- add_posterior_coefficients(add_seed(model, 8))
    expect_false(isTRUE(all.equal(third$posterior, first$posterior)))
  }
})

test_that("a seeded model draws the same under another kind of generator", {
  model <- add_seed(prepared_dfm()$object, 7)
  reference <- add_posterior_coefficients(model)$posterior

  # As on a cluster after parallel::clusterSetRNGStream().
  old <- RNGkind("L'Ecuyer-CMRG")
  set.seed(3)
  other <- add_posterior_coefficients(model)$posterior
  kind_after <- RNGkind()[1]
  RNGkind(old[1], old[2], old[3])

  expect_identical(other, reference)
  expect_identical(kind_after, "L'Ecuyer-CMRG")
})

test_that("a model without a seed draws from R's generator as it stands", {
  model <- prepared_dfm()$object
  model$model$seed <- NULL
  run <- function(state) {
    set.seed(state)
    add_posterior_coefficients(model)$posterior$factors$coeffs
  }

  expect_identical(run(5), run(5))
  expect_false(isTRUE(all.equal(run(5), run(6))))
})

test_that("a posterior_function is called as it is, seed included", {
  model <- add_seed(prepared_dfm()$object, 9)
  seen <- NULL
  add_posterior_coefficients(model, posterior_function = function(x) {
    seen <<- x$model$seed
    x
  })
  expect_identical(seen, 9L)
})

test_that("add_initial_values() gives every model of a list a seed of its own", {
  models <- create_dfmodel(x = sim_dfm()$x, p = 1:2, n = 1:2,
                           iterations = 20, burnin = 10)
  models <- add_initial_values(add_priors(models))

  seeds <- vapply(models, function(m) m$model$seed, integer(1))
  expect_length(unique(seeds), length(models))
})
