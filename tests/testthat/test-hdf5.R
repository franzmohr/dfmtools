# Writing a factor model into a BayesTS model file and reading it back.
#
# Three things are checked. A model comes back as it was written, for every
# sampler. A model read back draws what the original draws, which is what the
# translation is for: the file is in the sampler's names and ordering, and a
# mistake in either direction would change a prior without failing. And the
# file is in that ordering, checked against a permutation worked out by hand
# rather than against the code that computes it. The tree writer itself is
# bvartools'.

skip_if_not_installed("hdf5r")

# Lists in the order of their names, all the way down, since a file returns its
# groups alphabetically.
sorted_list <- function(x) {
  if (!is.list(x) || is.data.frame(x)) {
    return(x)
  }
  result <- lapply(x[order(names(x))], sorted_list)
  attributes(result) <- c(attributes(result)["names"],
                          attributes(x)[setdiff(names(attributes(x)), "names")])
  result
}

expect_same_model <- function(back, object) {
  # The writer gives a prior on the loadings without a mean the zero mean the
  # samplers imply, since BayesTS reads a prior only where it finds one.
  lambda <- object[["priors"]][["lambda"]]
  if (!is.null(lambda) && is.null(lambda[["mu"]])) {
    object[["priors"]][["lambda"]][["mu"]] <- matrix(0, nrow(lambda[["vinv"]]))
  }
  expect_identical(class(back), class(object))
  for (i in c("model", "data", "priors", "initial", "posterior")) {
    expect_equal(sorted_list(back[[i]]), sorted_list(object[[i]]), info = i)
  }
}

round_trip <- function(object) {
  path <- tempfile(fileext = ".h5")
  write_to_hdf5(object, filename = path)
  bvartools::read_model_from_hdf5(path)
}

# Two factors, so that the loadings' two orderings differ.
fitted_models <- function() {
  cached <- getOption("dfmtools.test.hdf5_fitted")
  if (!is.null(cached)) return(cached)
  prepared <- list(
    "DfmNormalGamma" = prepared_dfm(n = 2)$object,
    "DfmNormalStochvol" = prepared_dfm_sv(n = 2)$object,
    "DfmTvpGamma" = prepared_dfm_tvp(n = 2)$object,
    "DfmTvpStochvol" = prepared_dfm_tvp_sv(n = 2)$object
  )
  fitted <- lapply(prepared, function(object) {
    object <- add_posterior_coefficients(object)
    add_posterior_loglik(object)
  })
  result <- list("prepared" = prepared, "fitted" = fitted)
  options(dfmtools.test.hdf5_fitted = result)
  result
}

test_that("the permutation of the free loadings is the sampler's", {
  # Four series, two factors. R runs down the columns of lower.tri(): (2,1),
  # (3,1), (4,1), then (3,2), (4,2). The sampler runs along the rows: (2,1),
  # then (3,1), (3,2), then (4,1), (4,2) -- which are R's 1, 2, 4, 3, 5.
  expect_identical(.lambda_row_major_order(4L, 2L), c(1L, 2L, 4L, 3L, 5L))
  expect_identical(.lambda_row_major_order(4L, 1L), 1:3)
})

test_that("every factor model comes back as it was written", {
  models <- fitted_models()
  for (algorithm in names(models[["fitted"]])) {
    for (object in list(models[["prepared"]][[algorithm]],
                        models[["fitted"]][[algorithm]])) {
      expect_same_model(round_trip(object), object)
    }
  }
})

test_that("a model read back draws what the original draws", {
  models <- fitted_models()[["prepared"]]
  for (algorithm in names(models)) {
    object <- models[[algorithm]]
    # A prior that reads differently in the two orderings: a diagonal does not.
    n_lambda <- nrow(object[["priors"]][["lambda"]][["vinv"]])
    object[["priors"]][["lambda"]][["vinv"]] <- diag(seq_len(n_lambda) / 10, n_lambda)
    object[["priors"]][["lambda"]][["mu"]] <- matrix(seq_len(n_lambda) / 5)

    original <- add_posterior_coefficients(object)
    restored <- add_posterior_coefficients(round_trip(object))

    expect_equal(unclass(restored[["posterior"]][["lambda"]][["coeffs"]]),
                 unclass(original[["posterior"]][["lambda"]][["coeffs"]]),
                 info = algorithm)
    expect_equal(unclass(restored[["posterior"]][["factors"]][["coeffs"]]),
                 unclass(original[["posterior"]][["factors"]][["coeffs"]]),
                 info = algorithm)
  }
})

test_that("the file holds the sampler's names and ordering", {
  object <- fitted_models()[["prepared"]][["DfmTvpGamma"]]
  object[["priors"]][["lambda"]][["mu"]] <- matrix(c(1, 2, 3, 4, 5))
  object[["priors"]][["lambda"]][["shape"]] <- c(1, 2, 3, 4, 5)
  path <- tempfile(fileext = ".h5")
  write_to_hdf5(object, filename = path)

  tree <- bvartools::read_bayests_tree(path)
  spec <- tree[["model"]][[".attributes"]]

  expect_identical(spec[["k"]], 4L)
  expect_identical(spec[["n_factors"]], 2L)
  # BayesTS would read these as a VAR's exogenous variables and deterministic
  # terms.
  expect_false(any(c("m", "n") %in% names(spec)))

  expect_equal(as.vector(tree[["priors"]][["lambda"]][["mu"]]), c(1, 2, 4, 3, 5))
  expect_equal(as.vector(tree[["priors"]][["lambda"]][["shape"]]), c(1, 2, 4, 3, 5))
  expect_true(all(c("u_sigma", "v_sigma") %in% names(tree[["priors"]])))
  expect_true("v_inv" %in% names(tree[["priors"]][["lambda"]]))
  expect_equal(unclass(tree[["data"]][["train"]][["y"]]), unclass(object[["data"]][["x"]]),
               ignore_attr = TRUE)

  initial <- object[["initial"]]
  expect_equal(tree[["initial"]][["lambda"]], initial[["lambda"]][c(1, 2, 4, 3, 5), ],
               ignore_attr = TRUE)
  expect_equal(as.vector(tree[["initial"]][["u_sigma_inv"]]), diag(initial[["uinv"]]))
})

test_that("a FAVAR comes back as it was written, its identification included", {
  sim <- make_favar_sample(tt = 80)
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = 1, n = 1,
                             iterations = 60, burnin = 20)
  model <- add_priors(model)
  model <- add_initial_values(model)
  expect_same_model(round_trip(model), model)

  model <- add_posterior_coefficients(model)
  state <- c(colnames(model[["data"]][["x"]])[1], colnames(model[["data"]][["y"]]))
  restrictions <- data.frame(impulse = state[1], response = state[1], sign = 1)
  set.seed(3)
  identified <- add_sign_zero_restrictions(model, restrictions, max_tries = 20)

  back <- round_trip(identified)
  expect_same_model(back, identified)
  expect_s3_class(back[["model"]][["sign_zero_restrictions"]][["restrictions"]],
                  "data.frame")

  tree <- bvartools::read_bayests_tree({
    path <- tempfile(fileext = ".h5")
    write_to_hdf5(model, filename = path)
    path
  })
  expect_identical(tree[["model"]][[".attributes"]][["n_obs_factors"]],
                   model[["model"]][["n_obs"]])
  expect_equal(unclass(tree[["data"]][["train"]][["f_obs"]]),
               unclass(model[["data"]][["y"]]), ignore_attr = TRUE)
})

test_that("a file of the command line, without an R class, is read as a dfmodel", {
  object <- fitted_models()[["fitted"]][["DfmNormalGamma"]]
  tree <- .factor_model_tree(object)
  tree[["model"]][[".attributes"]][["rclass"]] <- NULL
  path <- tempfile(fileext = ".h5")
  bvartools::write_bayests_tree(tree, filename = path)

  back <- bvartools::read_model_from_hdf5(path)

  expect_s3_class(back, "dfmodel")
  expect_equal(back[["posterior"]][["lambda"]][["coeffs"]],
               object[["posterior"]][["lambda"]][["coeffs"]])
})

test_that("factor models and a VAR go through a folder in one modellist", {
  data("e1", package = "bvartools")
  var <- bvartools::create_bvarmodel(diff(log(e1)) * 100, p = 1,
                                     deterministic = "const",
                                     iterations = 20, burnin = 10)
  var <- bvartools::add_priors(var, coef = list(v_i = 1),
                               sigma = list(df = 3, scale = 1))
  var <- bvartools::add_initial_values(var)

  dfms <- create_dfmodel(x = sim_dfm()$x, p = 1:2, n = 1,
                         iterations = 20, burnin = 10)
  dfms <- add_initial_values(add_priors(dfms))
  models <- structure(list(dfms[[1]], var, dfms[[2]]), class = c("modellist", "list"))

  folder <- tempfile()
  dir.create(folder)
  write_to_hdf5(models, folder = folder)
  back <- bvartools::read_models_from_folder(folder)

  expect_s3_class(back, "modellist")
  expect_identical(unname(vapply(back, function(x) class(x)[1], "")),
                   c("dfmodel", "bvarmodel", "dfmodel"))
  expect_identical(back[[1]][["model"]][["p"]], 1L)
  expect_identical(back[[3]][["model"]][["p"]], 2L)
  expect_same_model(back[[3]], dfms[[2]])
})
