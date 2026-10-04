#' Add Forecast Errors to a Factor Model
#'
#' Calculates the forecast errors of a dynamic factor model or a factor augmented VAR, the
#' counterpart of \code{\link[bvartools]{add_forecast_errors.bvarmodel}} for a VAR.
#'
#' @param object an object of class \code{'dfmodel'} or \code{'favarmodel'}, usually the result
#' of a call to \code{\link{add_posterior_forecasts.dfmodel}} or
#' \code{\link{add_posterior_forecasts.favarmodel}}.
#' @param test_sample a time-series object used as test data. If \code{NULL} (default), the values
#' the model carries in \code{data$test} are used. For a dynamic factor model it holds the observed
#' series, one column per column of \code{data$x}. For a factor augmented VAR it holds the panel
#' series followed by the observed variables, the order of its forecast; columns named after all
#' of them are put into that order.
#' @param ... arguments passed forward to method.
#'
#' @details The realised values are matched to the forecast periods by time, the first forecast
#' period being the one after the estimation sample, and only the unbroken run of complete periods
#' from there is used. A test sample that does not reach the first forecast period leaves the model
#' unchanged.
#'
#' The forecasts are on the scale the model was estimated on, so the realised panel series are
#' put on it too -- centred and scaled with the moments of the estimation sample where
#' \code{normalize_x = TRUE} -- and the errors are on that scale, as the predictive likelihood of
#' \code{\link{add_predictive_loglik.dfmodel}} is. The observed variables of a factor augmented VAR
#' are never normalised.
#'
#' @return The object in \code{object} with \code{posterior$forecast$errors} added, a
#' \code{\link[coda]{mcmc}} object with one row per draw and one column per variable and realised
#' period, in the order of \code{posterior$forecast$forecasts}, which it is stored beside. The
#' realised values are stored in \code{data$test}, so the model can be scored again without them.
#'
#' @examples
#'
#' data("bem_dfmdata")
#' train <- window(bem_dfmdata, end = c(2014, 4))
#'
#' model <- create_dfmodel(x = train, p = 1, n = 1,
#'                         iterations = 20, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#' model <- add_forecast_input(model, n_ahead = 4)
#' model <- add_posterior_forecasts(model)
#'
#' model <- add_forecast_errors(model, test_sample = bem_dfmdata)
#' dim(get_forecast_errors(model))
#'
#' @export
#' @rdname add_forecast_errors.dfmodel
add_forecast_errors.dfmodel <- function(object, test_sample = NULL, ...) {

  realised <- if (is.null(test_sample)) {
    .realised_series(object)
  } else {
    .align_test_series(object, test_sample)
  }
  if (is.null(realised)) {
    return(object)
  }

  object[["data"]][["test"]][["x"]] <- realised
  .with_forecast_errors(object, realised)
}

#' @export
#' @rdname add_forecast_errors.dfmodel
add_forecast_errors.favarmodel <- function(object, test_sample = NULL, ...) {

  m <- object[["model"]][["m"]]
  n_obs <- object[["model"]][["n_obs"]]

  if (is.null(test_sample)) {
    x <- object[["data"]][["test"]][["x"]]
    y <- object[["data"]][["test"]][["y"]]
    if (is.null(x) || is.null(y)) {
      stop("No test sample was given and the object carries none in data$test. Pass ",
           "'test_sample'.", call. = FALSE)
    }
    realised <- cbind(as.matrix(x), as.matrix(y))
  } else {
    names <- c(colnames(object[["data"]][["x"]]), colnames(object[["data"]][["y"]]))
    if (!is.null(colnames(test_sample)) && all(names %in% colnames(test_sample))) {
      test_sample <- test_sample[, names, drop = FALSE]
    }
    if (NCOL(test_sample) != m + n_obs) {
      stop("Argument 'test_sample' has ", NCOL(test_sample), " columns, but the forecast ",
           "of the model has ", m + n_obs, ": the ", m, " panel series followed by the ",
           n_obs, " observed variables.", call. = FALSE)
    }
    realised <- .align_test_series(object, test_sample, width = m + n_obs, normalise = FALSE)
    if (is.null(realised)) {
      return(object)
    }
    realised[, seq_len(m)] <- .normalise_like_train(object, realised[, seq_len(m), drop = FALSE])
  }

  object[["data"]][["test"]][["x"]] <- realised[, seq_len(m), drop = FALSE]
  object[["data"]][["test"]][["y"]] <- realised[, m + seq_len(n_obs), drop = FALSE]
  .with_forecast_errors(object, realised)
}


# The model with posterior$forecast$errors: what was realised less each draw's
# forecast of it, over the periods `realised` has rows for.
.with_forecast_errors <- function(object, realised) {

  forecasts <- object[["posterior"]][["forecast"]][["forecasts"]]
  if (is.null(forecasts)) {
    stop("Object does not contain forecasts. Use add_posterior_forecasts first.", call. = FALSE)
  }

  realised <- unclass(as.matrix(realised))
  width <- ncol(realised)
  h <- nrow(realised)
  draws <- nrow(forecasts)
  mc_stats <- coda::mcpar(forecasts)

  errors <- t(matrix(t(realised), h * width, draws)) -
    unclass(as.matrix(forecasts))[, seq_len(h * width), drop = FALSE]
  object[["posterior"]][["forecast"]][["errors"]] <-
    coda::mcmc(errors, start = mc_stats[1], end = mc_stats[2], thin = mc_stats[3])

  return(object)
}


#' Get Forecast Errors of a Factor Model
#'
#' Returns the forecast errors \code{\link{add_forecast_errors.dfmodel}} added to a model.
#'
#' @param object an object of class \code{'dfmodel'} or \code{'favarmodel'}.
#' @param ... not used.
#'
#' @return The \code{\link[coda]{mcmc}} object in \code{posterior$forecast$errors}, or
#' \code{NULL} where there is none.
#'
#' @export
#' @rdname get_forecast_errors.dfmodel
get_forecast_errors.dfmodel <- function(object, ...) {
  return(object[["posterior"]][["forecast"]][["errors"]])
}

#' @export
#' @rdname get_forecast_errors.dfmodel
get_forecast_errors.favarmodel <- function(object, ...) {
  return(object[["posterior"]][["forecast"]][["errors"]])
}
