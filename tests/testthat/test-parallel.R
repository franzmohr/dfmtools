# Lists of models on several cores. The list methods and their argument 'cores'
# are bvartools'; what is tested here is that they work for the models of this
# package, whose methods a worker only has once it has loaded dfmtools.
#
# The workers load the installed dfmtools, so these tests say something about
# the code under test only when that is what is installed -- under R CMD check,
# not under devtools::load_all().
skip_if_workers_load_other_version <- function() {
  skip_on_cran()
  skip_if(isNamespaceLoaded("pkgload") && pkgload::is_dev_package("dfmtools"),
          "workers would load the installed dfmtools, not the one under test")
}

# Two lag orders each: a 'modellist' of two models.
dfm_list <- function() {
  models <- create_dfmodel(x = sim_dfm()$x, p = 1:2, n = 1,
                           iterations = 20, burnin = 10)
  models <- add_priors(models)
  set.seed(41)
  add_initial_values(models)
}

favar_list <- function() {
  sample <- make_favar_sample(tt = 60, n_x = 5)
  models <- create_favarmodel(x = sample$x, y = sample$yts, p = 1:2, n = 1,
                              normalize_x = FALSE, iterations = 20, burnin = 10)
  models <- add_priors(models)
  set.seed(42)
  add_initial_values(models)
}

posteriors <- function(object, element = NULL) {
  lapply(object, function(model) {
    if (is.null(element)) model$posterior else model$posterior[[element]]
  })
}

# Evaluates 'fun(object, ...)' one model after the other in a separate process
# that runs with one thread, as the workers do. The thread count changes the
# last digits of some results -- those of the FAVAR sampler among them -- and a
# chain carries the difference forward, so this session is no reference.
on_one_thread <- function(fun, object, ...) {
  vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
            "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  old <- Sys.getenv(vars, unset = NA, names = TRUE)
  do.call(Sys.setenv, stats::setNames(as.list(rep("1", length(vars))), vars))
  cl <- tryCatch(parallel::makeCluster(1), finally = {
    for (v in vars) {
      if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
    }
  })
  on.exit(parallel::stopCluster(cl), add = TRUE)

  attach_dfmtools <- function(paths) {
    .libPaths(paths)
    suppressPackageStartupMessages(library("dfmtools", character.only = TRUE))
    invisible(NULL)
  }
  environment(attach_dfmtools) <- baseenv()
  parallel::clusterCall(cl, attach_dfmtools, .libPaths())
  parallel::clusterCall(cl, fun, object, ...)[[1]]
}

no_errors <- function(object) {
  !any(vapply(object, function(model) isTRUE(model$error), logical(1)))
}

test_that("a list of dynamic factor models is simulated on two workers", {
  skip_if_workers_load_other_version()
  models <- dfm_list()

  parallel <- add_posterior_coefficients(models, cores = 2)
  expect_s3_class(parallel, "modellist")
  expect_true(no_errors(parallel))
  expect_identical(posteriors(parallel),
                   posteriors(on_one_thread(add_posterior_coefficients, models)))

  loglik <- add_posterior_loglik(parallel, cores = 2)
  expect_identical(posteriors(loglik, "loglik"),
                   posteriors(on_one_thread(add_posterior_loglik, parallel), "loglik"))

  # Forecasts are not seeded per model; they draw from streams of the workers
  # that set.seed() fixes.
  set.seed(5)
  first <- add_posterior_forecasts(parallel, n_ahead = 2, cores = 2)
  set.seed(5)
  second <- add_posterior_forecasts(parallel, n_ahead = 2, cores = 2)
  expect_identical(posteriors(first, "forecast"), posteriors(second, "forecast"))
  expect_identical(lapply(posteriors(first, "forecast"), dim),
                   lapply(posteriors(add_posterior_forecasts(parallel, n_ahead = 2),
                                     "forecast"), dim))
})

test_that("a list of factor augmented VARs is simulated on two workers", {
  skip_if_workers_load_other_version()
  models <- favar_list()

  parallel <- add_posterior_coefficients(models, cores = 2)
  expect_s3_class(parallel, "modellist")
  expect_true(no_errors(parallel))
  expect_identical(posteriors(parallel),
                   posteriors(on_one_thread(add_posterior_coefficients, models)))

  loglik <- add_posterior_loglik(parallel, cores = 2)
  expect_true(no_errors(loglik))
  expect_identical(posteriors(loglik, "loglik"),
                   posteriors(on_one_thread(add_posterior_loglik, parallel), "loglik"))
})

test_that("a model of a list without a seed is given one before it is sent off", {
  skip_if_workers_load_other_version()
  models <- dfm_list()
  models[[2]]$model$seed <- NULL

  set.seed(8)
  first <- add_posterior_coefficients(models, cores = 2)
  set.seed(8)
  second <- add_posterior_coefficients(models, cores = 2)

  expect_type(first[[2]]$model$seed, "integer")
  expect_identical(first[[2]]$posterior, second[[2]]$posterior)
})
