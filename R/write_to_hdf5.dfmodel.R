#' Export to and Import from HDF5 Files
#'
#' Writes a dynamic factor model or a factor augmented VAR into an HDF5 model
#' file that the BayesTS command line can run, and reads it back.
#'
#' @param object an object of class \code{'dfmodel'} or \code{'favarmodel'}.
#' @param tree the tree of a model file, as
#' \code{\link[bvartools]{read_bayests_tree}} returns it.
#' @param filename path to the file, in which output should be stored.
#' @param group the group the model's tree should hang under inside its file.
#' Defaults to \code{""}, the root of the file, which is where a file holding a
#' single model puts it.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The file is in the layout BayesTS reads, so a model written here can be
#' estimated, forecast and scored by the command line -- with
#' \code{\link[bvartools]{bayests_files}}, say -- and read back with
#' \code{\link[bvartools]{read_model_from_hdf5}}, which returns the
#' \code{'dfmodel'} or \code{'favarmodel'} it was written from. That layout is
#' not the one of the R object, and the writer translates:
#' \itemize{
#'   \item the number of observed series, \code{model$m}, is \code{k} in the
#'   file, the number of factors, \code{model$n}, is \code{n_factors}, and a
#'   FAVAR's observed variables, \code{model$n_obs}, are \code{n_obs_factors}.
#'   \strong{BayesTS reads \code{m} and \code{n} as the numbers of exogenous
#'   variables and deterministic terms of a VAR}, so they are not written under
#'   their R names;
#'   \item the panel \code{data$x} is \code{/data/train/y}, what the horizon
#'   realised, \code{data$test$x}, is \code{/data/test/y}, and the observed
#'   variables of a FAVAR, \code{data$y}, are \code{/data/train/f_obs};
#'   \item the deterministic terms, \code{data$deterministic}, are
#'   \code{/data/train/x}, and over the horizon \code{data$forecast$x} is
#'   \code{/data/forecast/x}. \strong{Their number is \code{/model/n}}, which
#'   BayesTS checks against the columns of \code{/data/train/x};
#'   \item the priors \code{u} and \code{v} are \code{/priors/u_sigma} and
#'   \code{/priors/v_sigma}, the starting precisions \code{uinv} and
#'   \code{vinv} are \code{/initial/u_sigma_inv} and \code{/initial/v_sigma_inv}
#'   -- stored as their diagonals where they are diagonal -- and a prior
#'   precision \code{vinv} is \code{v_inv};
#'   \item \strong{the free loadings of a dynamic factor model are put in the
#'   order of the sampler}, row by row, where R keeps them column by column.
#'   This applies to everything that has one entry per free loading: the prior,
#'   the starting values and, for time-varying loadings, the draws of their
#'   random walk's variance. The permutation is the one the samplers of this
#'   package use, so a model gives the same draws whether it is estimated here
#'   or from its file. A FAVAR's loadings are flattened row by row, as its
#'   sampler takes them.
#' }
#' A prior on the loadings without a mean is written with the zero mean the
#' samplers imply, since BayesTS reads a prior only where it finds one.
#'
#' Everything else -- the remaining specification, the draws, the forecasts and
#' the log-likelihoods, and the result of
#' \code{\link{add_sign_zero_restrictions.favarmodel}} -- keeps its R name.
#'
#' \code{from_bayests_tree()} undoes the translation. It is called by
#' \code{\link[bvartools]{read_model_from_hdf5}}, which loads this package for
#' a factor model's file whether or not it is attached, and is not meant to be
#' called directly.
#'
#' A list of models, factor models among them, is written with the
#' \code{'modellist'} method of \pkg{bvartools}, one file per model, and read
#' back with \code{\link[bvartools]{read_models_from_folder}}.
#'
#' @return \code{write_to_hdf5()} returns the path to the written file,
#' invisibly. \code{from_bayests_tree()} returns an object of class
#' \code{'dfmodel'} or \code{'favarmodel'}.
#'
#' @examples
#' data("bem_dfmdata")
#'
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1,
#'                         iterations = 20, burnin = 10)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#'
#' path <- tempfile(fileext = ".h5")
#' write_to_hdf5(model, filename = path)
#' model <- bvartools::read_model_from_hdf5(path)
#'
#' @export
write_to_hdf5.dfmodel <- function(object, filename, group = "", ...) {
  bvartools::write_bayests_tree(.factor_model_tree(object), filename = filename,
                                group = group)
}

#' @rdname write_to_hdf5.dfmodel
#' @export
write_to_hdf5.favarmodel <- function(object, filename, group = "", ...) {
  bvartools::write_bayests_tree(.factor_model_tree(object), filename = filename,
                                group = group)
}

#' @rdname write_to_hdf5.dfmodel
#' @export
from_bayests_tree.dfmodel <- function(tree, ...) {
  .factor_model_from_tree(tree)
}

#' @rdname write_to_hdf5.dfmodel
#' @export
from_bayests_tree.favarmodel <- function(tree, ...) {
  .factor_model_from_tree(tree)
}


# The names a factor model's specification has in R and in a model file. Each
# R name is written only under the file's name, which is the point: m and n
# mean something else to BayesTS.
.factor_spec_names <- c("m" = "k", "n" = "n_factors", "n_obs" = "n_obs_factors")

# The prior groups of a normal block, whose precision R calls vinv and BayesTS
# v_inv: the loadings, the transition, and the coefficients of the
# deterministic terms in the measurement equation and, for a FAVAR, in the
# observed variables'.
.normal_blocks <- c("lambda", "a", "c", "c_obs")

# The elements of the specification BayesTS reads as integers.
.factor_spec_counts <- c("k", "n_factors", "n_obs_factors", "p", "iterations",
                         "burnin", "thin", "h", "seed")

# The panel's normalisation, the centres and scales scale() attaches to it,
# which a model needs to be read on the scale of its data. They are kept beside
# the panel rather than on it, as /data/scaling, which no sampler reads.
.scaling_names <- c("center" = "scaled:center", "scale" = "scaled:scale")

# What a factor model is translated with: whether it is a FAVAR, and the
# permutation of its free loadings, which a FAVAR does not need.
.factor_translation <- function(model) {
  favar <- identical(model[["algorithm"]], "FavarNormalWishart")
  list("favar" = favar,
       "tvp" = isTRUE(model[["tvp"]]),
       "order" = if (favar) NULL else .lambda_row_major_order(model[["m"]], model[["n"]]))
}

# A value with one entry per free loading, put on the sampler's ordering, or
# back on R's. A value of another length -- a scalar the sampler recycles, or a
# wrong size the sampler should refuse by its own size -- is left as it is,
# which is what the bindings do. Dimensions are kept, so a column stays a
# column.
.to_core_order <- function(x, order) {
  if (is.null(x) || is.null(order) || length(x) != length(order)) {
    return(x)
  }
  x[] <- as.vector(x)[order]
  x
}

.to_r_order <- function(x, order) {
  if (is.null(x) || is.null(order) || length(x) != length(order)) {
    return(x)
  }
  value <- as.vector(x)
  value[order] <- as.vector(x)
  x[] <- value
  x
}

# The same for a matrix with one row (`rows`), one column (`columns`) or both
# per free loading.
.permute_loadings <- function(x, order, rows = FALSE, columns = FALSE, back = FALSE) {
  if (is.null(x) || is.null(order) || is.null(dim(x))) {
    return(x)
  }
  index <- if (back) order(order) else order
  if (rows && nrow(x) == length(order)) {
    x <- x[index, , drop = FALSE]
  }
  if (columns && ncol(x) == length(order)) {
    # Subsetting an mcmc object by columns keeps its mcpar; a plain matrix keeps
    # what it has.
    x <- x[, index, drop = FALSE]
  }
  x
}

# A diagonal precision as the vector of its diagonal, a column as the file
# keeps it, and back.
.diagonal_column <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  matrix(diag(as.matrix(x)), ncol = 1)
}

.diagonal_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  diag(as.vector(x), nrow = length(x))
}

# A list that may hold a data frame -- the restrictions of
# add_sign_zero_restrictions() -- as a tree, and back. A data frame is a group
# of its columns, which HDF5 lists in alphabetical order, so the order is kept
# as an attribute.
.to_tree_value <- function(x) {
  if (is.data.frame(x)) {
    return(c(lapply(as.list(x), .to_tree_value),
             list(".attributes" = list("rclass" = "data.frame",
                                       "columns" = names(x)))))
  }
  if (is.list(x)) {
    return(lapply(x, .to_tree_value))
  }
  x
}

.from_tree_value <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  attrs <- x[[".attributes"]]
  x[[".attributes"]] <- NULL
  x <- lapply(x, .from_tree_value)
  if (identical(attrs[["rclass"]], "data.frame")) {
    return(as.data.frame(x[attrs[["columns"]]], stringsAsFactors = FALSE,
                         optional = TRUE))
  }
  x
}

# A factor model as the tree of its model file.
.factor_model_tree <- function(object) {

  model <- object[["model"]]
  if (is.null(model[["algorithm"]])) {
    stop("Element 'model$algorithm' is missing. Was the object produced by ",
         "create_dfmodel() or create_favarmodel()?")
  }
  tr <- .factor_translation(model)
  order <- tr[["order"]]

  # Specification ----
  attributes <- list()
  nested <- list()
  for (i in names(model)) {
    value <- model[[i]]
    if (is.null(value)) {
      next
    }
    name <- if (i %in% names(.factor_spec_names)) .factor_spec_names[[i]] else i
    if (is.list(value)) {
      nested[[name]] <- .to_tree_value(value)
    } else {
      # BayesTS reads the counts as integers, and R keeps some of them as
      # doubles -- the number of factors of a model of a grid, for one.
      if (name %in% .factor_spec_counts && is.numeric(value) &&
          all(value == round(value))) {
        value <- as.integer(value)
      }
      attributes[[name]] <- value
    }
  }
  # The number of deterministic terms, which is what BayesTS calls n. The core
  # takes it from the columns of the terms when it is called from R, but reads
  # it from /model/n in a file, and refuses /data/train/x when they disagree.
  if (!is.null(object[["data"]][["deterministic"]])) {
    attributes[["n"]] <- ncol(object[["data"]][["deterministic"]])
  }
  attributes[["rclass"]] <- class(object)
  # The package bvartools loads to read the file back with the methods below.
  attributes[["rpackage"]] <- "dfmtools"
  tree_model <- c(nested, list(".attributes" = attributes))

  # Data ----
  data <- object[["data"]]
  scaling <- lapply(.scaling_names, function(i) attr(data[["x"]], i))
  # The deterministic terms are the regressors of the measurement equation,
  # /data/train/x, and over the horizon /data/forecast/x. What the horizon
  # realised is /data/test/y for the panel and, for a FAVAR, /data/test/f_obs
  # for its observed variables, which BayesTS does not read.
  tree_data <- list(
    "train" = list("y" = data[["x"]], "f_obs" = data[["y"]],
                   "x" = data[["deterministic"]]),
    "forecast" = data[["forecast"]],
    "scaling" = if (!all(vapply(scaling, is.null, logical(1)))) scaling,
    "test" = if (!is.null(data[["test"]][["x"]]) || !is.null(data[["test"]][["y"]])) {
      list("y" = data[["test"]][["x"]], "f_obs" = data[["test"]][["y"]])
    }
  )

  # Priors ----
  priors <- object[["priors"]]
  tree_priors <- NULL
  if (length(priors) > 0) {
    tree_priors <- list()
    for (i in names(priors)) {
      block <- priors[[i]]
      name <- switch(i, "u" = "u_sigma", "v" = "v_sigma", i)
      if (i %in% .normal_blocks) {
        names(block)[names(block) == "vinv"] <- "v_inv"
      }
      # BayesTS reads a block's prior only where it finds a mean. The bindings
      # imply a zero one.
      if (i %in% c("c", "c_obs") && is.null(block[["mu"]]) && !is.null(block[["v_inv"]])) {
        block[["mu"]] <- matrix(0, nrow(block[["v_inv"]]))
      }
      if (i == "lambda" && !tr[["favar"]]) {
        if (is.null(block[["mu"]]) && !is.null(order)) {
          block[["mu"]] <- matrix(0, length(order))
        }
        block[["v_inv"]] <- .permute_loadings(block[["v_inv"]], order,
                                              rows = TRUE, columns = TRUE)
        for (j in c("mu", "shape", "rate", "omega_v")) {
          block[[j]] <- .to_core_order(block[[j]], order)
        }
      }
      tree_priors[[name]] <- .to_tree_value(block)
    }
  }

  # Initial values ----
  initial <- object[["initial"]]
  tree_initial <- NULL
  if (length(initial) > 0) {
    tree_initial <- initial
    tree_initial[["uinv"]] <- NULL
    tree_initial[["vinv"]] <- NULL
    tree_initial[["u_sigma_inv"]] <- .diagonal_column(initial[["uinv"]])
    # A FAVAR's state innovation precision is a matrix, not a diagonal.
    tree_initial[["v_sigma_inv"]] <- if (tr[["favar"]]) {
      initial[["vinv"]]
    } else {
      .diagonal_column(initial[["vinv"]])
    }
    if (tr[["favar"]]) {
      # (k - n) x n_state in R, row by row in the sampler.
      if (!is.null(initial[["lambda"]])) {
        tree_initial[["lambda"]] <- matrix(as.vector(t(initial[["lambda"]])), ncol = 1)
      }
    } else if (tr[["tvp"]]) {
      tree_initial[["lambda"]] <- .permute_loadings(initial[["lambda"]], order, rows = TRUE)
      tree_initial[["lambda_sigma_inv"]] <-
        .permute_loadings(initial[["lambda_sigma_inv"]], order, rows = TRUE, columns = TRUE)
      tree_initial[["lambda_init"]] <- .to_core_order(initial[["lambda_init"]], order)
    } else {
      tree_initial[["lambda"]] <- .to_core_order(initial[["lambda"]], order)
    }
  }

  # Posterior ----
  # Under the sampler's names already. What has one entry per free loading --
  # the variance of a time-varying loading's random walk and, under the
  # non-centred prior, its signed standard deviation and the ordinate at zero
  # -- is a column per loading in R and goes on the sampler's order.
  posterior <- object[["posterior"]]
  if (!is.null(posterior[["lambda"]]) && !is.null(order)) {
    for (j in c("sigma", "omega", "omega_log_zero")) {
      posterior[["lambda"]][[j]] <- .permute_loadings(posterior[["lambda"]][[j]],
                                                      order, columns = TRUE)
    }
  }

  list("model" = tree_model,
       "data" = tree_data,
       "priors" = tree_priors,
       "initial" = tree_initial,
       "posterior" = posterior)
}

# The tree of a factor model's file as the model.
.factor_model_from_tree <- function(tree) {

  attributes <- tree[["model"]][[".attributes"]]
  model_class <- attributes[["rclass"]]
  if (is.null(model_class)) {
    model_class <- class(tree)
  }
  attributes[["rclass"]] <- NULL
  attributes[["rpackage"]] <- NULL

  # Specification ----
  model <- list()
  # m and n of the file are BayesTS's: no exogenous variables, and the number
  # of deterministic terms, which the terms themselves say. The model's own
  # come from k and n_factors.
  attributes[["m"]] <- NULL
  attributes[["n"]] <- NULL
  for (i in names(attributes)) {
    name <- names(.factor_spec_names)[.factor_spec_names == i]
    model[[if (length(name) == 1) name else i]] <- attributes[[i]]
  }
  # A dynamic factor model has no observed variables and no n_obs, where a file
  # of the command line may still say n_obs_factors = 0.
  if (!"favarmodel" %in% model_class) {
    model[["n_obs"]] <- NULL
  }
  nested <- tree[["model"]]
  nested[[".attributes"]] <- NULL
  for (i in names(nested)) {
    model[[i]] <- .from_tree_value(nested[[i]])
  }

  tr <- .factor_translation(model)
  order <- tr[["order"]]

  # Data ----
  train <- tree[["data"]][["train"]]
  data <- list("x" = train[["y"]])
  for (i in names(.scaling_names)) {
    value <- tree[["data"]][["scaling"]][[i]]
    if (!is.null(value)) {
      names(value) <- colnames(data[["x"]])
      attr(data[["x"]], .scaling_names[[i]]) <- value
    }
  }
  if (!is.null(train[["f_obs"]])) {
    data[["y"]] <- train[["f_obs"]]
  }
  if (!is.null(train[["x"]])) {
    data[["deterministic"]] <- train[["x"]]
  }
  if (!is.null(tree[["data"]][["forecast"]])) {
    data[["forecast"]] <- tree[["data"]][["forecast"]]
    # Named after the terms, as add_forecast_input() names them; a dataset
    # without series attributes keeps no names.
    if (!is.null(data[["forecast"]][["x"]]) &&
        ncol(data[["forecast"]][["x"]]) == length(model[["deterministic"]])) {
      colnames(data[["forecast"]][["x"]]) <- model[["deterministic"]]
    }
  }
  test <- tree[["data"]][["test"]]
  if (!is.null(test)) {
    data[["test"]] <- list("x" = test[["y"]], "y" = test[["f_obs"]])
    data[["test"]] <- data[["test"]][!vapply(data[["test"]], is.null, logical(1))]
  }

  result <- list("data" = data, "model" = model)

  # Priors ----
  if (!is.null(tree[["priors"]])) {
    priors <- list()
    for (name in names(tree[["priors"]])) {
      block <- .from_tree_value(tree[["priors"]][[name]])
      i <- switch(name, "u_sigma" = "u", "v_sigma" = "v", name)
      if (i == "lambda" && !tr[["favar"]]) {
        block[["v_inv"]] <- .permute_loadings(block[["v_inv"]], order, rows = TRUE,
                                              columns = TRUE, back = TRUE)
        for (j in c("mu", "shape", "rate", "omega_v")) {
          block[[j]] <- .to_r_order(block[[j]], order)
        }
      }
      if (i %in% .normal_blocks) {
        names(block)[names(block) == "v_inv"] <- "vinv"
      }
      priors[[i]] <- block
    }
    result[["priors"]] <- priors
  }

  # Initial values ----
  if (!is.null(tree[["initial"]])) {
    initial <- tree[["initial"]]
    initial[["u_sigma_inv"]] <- NULL
    initial[["v_sigma_inv"]] <- NULL
    initial[["uinv"]] <- .diagonal_matrix(tree[["initial"]][["u_sigma_inv"]])
    initial[["vinv"]] <- if (tr[["favar"]]) {
      tree[["initial"]][["v_sigma_inv"]]
    } else {
      .diagonal_matrix(tree[["initial"]][["v_sigma_inv"]])
    }
    if (tr[["favar"]]) {
      if (!is.null(initial[["lambda"]])) {
        initial[["lambda"]] <- matrix(as.vector(initial[["lambda"]]),
                                      nrow = model[["m"]] - model[["n"]],
                                      byrow = TRUE)
      }
    } else if (tr[["tvp"]]) {
      initial[["lambda"]] <- .permute_loadings(initial[["lambda"]], order,
                                               rows = TRUE, back = TRUE)
      initial[["lambda_sigma_inv"]] <- .permute_loadings(initial[["lambda_sigma_inv"]],
                                                         order, rows = TRUE,
                                                         columns = TRUE, back = TRUE)
      initial[["lambda_init"]] <- .to_r_order(initial[["lambda_init"]], order)
    } else {
      initial[["lambda"]] <- .to_r_order(initial[["lambda"]], order)
    }
    result[["initial"]] <- initial
  }

  # Posterior ----
  posterior <- tree[["posterior"]]
  if (!is.null(posterior[["lambda"]]) && !is.null(order)) {
    for (j in c("sigma", "omega", "omega_log_zero")) {
      posterior[["lambda"]][[j]] <- .permute_loadings(posterior[["lambda"]][[j]],
                                                      order, columns = TRUE, back = TRUE)
    }
  }
  # The coefficients of the deterministic terms, named "<series>.<term>" as
  # add_posterior_coefficients() names them. Draws keep no names in a file.
  for (i in c("c", "c_obs")) {
    draws <- posterior[[i]][["coeffs"]]
    series <- colnames(data[[if (i == "c") "x" else "y"]])
    if (!is.null(draws) && !is.null(model[["deterministic"]])) {
      labels <- .deterministic_draw_names(series, model[["deterministic"]])
      if (length(labels) == ncol(draws)) {
        colnames(posterior[[i]][["coeffs"]]) <- labels
      }
    }
  }
  if (!is.null(posterior)) {
    result[["posterior"]] <- posterior
  }

  class(result) <- model_class
  result
}
