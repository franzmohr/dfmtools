# What the factor augmented VAR's methods refuse, and what they do when the
# sampler itself fails.
#
# These are the error branches of the three methods that carry a model on from
# add_posterior_coefficients(). They were the least covered code in the package
# -- which is also where an audit found the prior arguments going unchecked, so
# the two agree about where the thin spot is.

favar_drawn <- function(iterations = 20, burnin = 10) {
  sim <- make_favar_sample(tt = 60, n_x = 5)
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1,
                             normalize_x = FALSE,
                             iterations = iterations, burnin = burnin)
  model <- add_initial_values(add_priors(model))
  add_posterior_coefficients(model)
}

test_that("an algorithm the methods do not know is refused by name", {

  drawn <- favar_drawn()
  wrong <- drawn
  wrong$model$algorithm <- "FavarSomethingElse"

  expect_error(add_posterior_forecasts(add_forecast_input(wrong, n_ahead = 2)), "not supported")
  expect_error(add_posterior_loglik(wrong), "not supported")

  # add_posterior_coefficients is the one step that carries on rather than
  # stopping, so the same mistake comes back as a flag on the model.
  refused <- suppressWarnings(add_posterior_coefficients(wrong))
  expect_true(refused$error)
  expect_null(refused$posterior)
})

test_that("a model with no algorithm says which function should have set one", {

  drawn <- favar_drawn()
  without <- drawn
  without$model$algorithm <- NULL

  expect_error(add_posterior_forecasts(add_forecast_input(without, n_ahead = 2)),
               "create_favarmodel")
  expect_error(add_posterior_loglik(without), "create_favarmodel")
  expect_true(suppressWarnings(add_posterior_coefficients(without))$error)
})

test_that("the derived steps need the draws they are derived from", {

  sim <- make_favar_sample(tt = 60, n_x = 5)
  undrawn <- add_initial_values(add_priors(
    create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1, normalize_x = FALSE,
                      iterations = 20, burnin = 10)))

  expect_error(add_posterior_loglik(undrawn), "does not contain posterior draws")
  expect_error(add_posterior_forecasts(add_forecast_input(undrawn, n_ahead = 2)),
               "does not contain posterior draws")

  # The state path is a member of the posterior rather than something the
  # parameters imply, so its absence is its own message.
  drawn <- favar_drawn()
  no_state <- drawn
  no_state$posterior$factors <- NULL
  expect_error(add_posterior_loglik(no_state), "posterior draws of the state")
  expect_error(add_posterior_forecasts(add_forecast_input(no_state, n_ahead = 2)),
               "posterior draws of the factors")
})

test_that("a failing sampler leaves the model with error = TRUE and no posterior", {

  # A starting value of the wrong width is refused by the core's validator,
  # which is the failure add_posterior_coefficients() is meant to survive.
  sim <- make_favar_sample(tt = 60, n_x = 5)
  broken <- add_initial_values(add_priors(
    create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1, normalize_x = FALSE,
                      iterations = 20, burnin = 10)))
  broken$initial$uinv <- diag(1, 3)

  out <- suppressWarnings(add_posterior_coefficients(broken))
  expect_true(out$error)
  expect_null(out$posterior)
  expect_s3_class(out, "favarmodel")
})

test_that("a posterior_function replaces the sampler and its failure is caught", {

  sim <- make_favar_sample(tt = 60, n_x = 5)
  model <- add_initial_values(add_priors(
    create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1, normalize_x = FALSE,
                      iterations = 20, burnin = 10)))

  own <- add_posterior_coefficients(model, posterior_function = function(object) {
    object$posterior <- list(note = "mine")
    object
  })
  expect_identical(own$posterior$note, "mine")
  expect_s3_class(own, "favarmodel")

  failing <- suppressWarnings(add_posterior_coefficients(
    model, posterior_function = function(object) stop("no draws")))
  expect_true(failing$error)
})
