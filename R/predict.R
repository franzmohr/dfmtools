#' Predict Method for Factor Models
#'
#' Collects the forecasts of a dynamic factor model or a factor augmented VAR, the counterpart of
#' \code{\link[bvartools]{predict.bvarmodel}} for a VAR.
#'
#' @param object an object of class \code{'dfmodel'} or \code{'favarmodel'}, usually the result
#' of a call to \code{\link{add_posterior_forecasts.dfmodel}} or
#' \code{\link{add_posterior_forecasts.favarmodel}}.
#' @param n_ahead number of steps ahead at which to predict. If \code{NULL} (default), every
#' period simulated by \code{\link{add_posterior_forecasts.dfmodel}}, i.e. the horizon given to
#' \code{\link{add_forecast_input.dfmodel}}.
#' @param rescale logical. If \code{TRUE} (default), the normalisation that
#' \code{normalize_x = TRUE} of \code{\link{create_dfmodel}} and \code{\link{create_favarmodel}}
#' applied to the panel is undone, so that the forecasts and the training sample are in the
#' units of the data. A model created with \code{normalize_x = FALSE} is in those units already.
#' @param ... additional arguments.
#'
#' @details The forecasts are draws from the posterior predictive distribution, one per posterior
#' draw, as \code{\link{add_posterior_forecasts.dfmodel}} simulated them. Those are on the scale the
#' model was estimated on; \code{rescale} is what puts them back.
#'
#' The forecast of a factor augmented VAR covers its observed variables as well as its panel, so
#' its result has a column for each of them, after the panel series. The observed variables are
#' never normalised and \code{rescale} leaves them alone.
#'
#' @return An object of class \code{'bvarprd'}, as \code{\link[bvartools]{predict.bvarmodel}}
#' returns: a list with element \code{fcst}, an \eqn{h \times K \times S} array of the \eqn{S}
#' forecast draws, whose first two dimensions are named after the forecast periods and the
#' variables, and element \code{y}, the variables of the training sample. bvartools' \code{plot}
#' method draws it.
#'
#' @examples
#'
#' data("bem_dfmdata")
#'
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1,
#'                         iterations = 20, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#' model <- add_forecast_input(model, n_ahead = 4)
#' model <- add_posterior_forecasts(model)
#'
#' pred <- predict(model)
#' dim(pred$fcst)
#'
#' @export
#' @rdname predict.dfmodel
predict.dfmodel <- function(object, n_ahead = NULL, rescale = TRUE, ...) {
  .factor_predict(object, n_ahead, rescale, extra = NULL)
}

#' @export
#' @rdname predict.dfmodel
predict.favarmodel <- function(object, n_ahead = NULL, rescale = TRUE, ...) {
  .factor_predict(object, n_ahead, rescale, extra = object[["data"]][["y"]])
}


# The forecasts of a factor model as bvartools' predict.bvarmodel() returns
# those of a VAR. `extra` is a FAVAR's observed variables, which its forecast
# carries after the panel in every period and which are never normalised.
.factor_predict <- function(object, n_ahead, rescale, extra) {

  h <- object[["model"]][["h"]]
  if (is.null(h)) {
    stop("Missing specification of h in object$model$h. You might want to use\n",
         "function 'add_forecast_input' and then 'add_posterior_forecasts'\n",
         "before using this function.", call. = FALSE)
  }
  forecasts <- object[["posterior"]][["forecast"]][["forecasts"]]
  if (is.null(forecasts)) {
    stop("Missing element object$posterior$forecast$forecasts. You might want to use\n",
         "'add_posterior_forecasts'\nbefore this function.", call. = FALSE)
  }
  if (!is.logical(rescale) || length(rescale) != 1 || is.na(rescale)) {
    stop("Argument 'rescale' must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(n_ahead)) {
    n_ahead <- h
  }
  .check_whole_number(n_ahead, "n_ahead", 1)
  if (h < n_ahead) {
    warning("Argument 'n_ahead' is larger than the value in object$model$h.\n",
            "Limiting the output to the latter.", call. = FALSE)
    n_ahead <- h
  }

  panel <- object[["data"]][["x"]]
  m <- NCOL(panel)
  center <- rep(0, m)
  scale <- rep(1, m)
  if (rescale) {
    if (!is.null(attr(panel, "scaled:center"))) {
      center <- as.numeric(attr(panel, "scaled:center"))
    }
    if (!is.null(attr(panel, "scaled:scale"))) {
      scale <- as.numeric(attr(panel, "scaled:scale"))
    }
  }

  tsp_train <- stats::tsp(panel)
  train <- sweep(sweep(unclass(as.matrix(panel)), 2, scale, "*"), 2, center, "+")
  varnames <- colnames(panel)
  if (!is.null(extra)) {
    train <- cbind(train, unclass(as.matrix(extra)))
    varnames <- c(varnames, colnames(extra))
    center <- c(center, rep(0, NCOL(extra)))
    scale <- c(scale, rep(1, NCOL(extra)))
  }
  attr(train, "scaled:center") <- NULL
  attr(train, "scaled:scale") <- NULL
  dimnames(train) <- list(NULL, varnames)
  train <- stats::ts(train, start = tsp_train[1], frequency = tsp_train[3])

  k <- length(varnames)
  forecasts <- unclass(as.matrix(forecasts))
  draws <- nrow(forecasts)
  prd_time <- tsp_train[2] + seq_len(n_ahead) / tsp_train[3]

  result <- array(NA_real_, dim = c(n_ahead, k, draws))
  dimnames(result) <- list(as.character(prd_time), varnames, NULL)
  for (i in seq_len(draws)) {
    result[, , i] <- t(matrix(forecasts[i, ], k)[, seq_len(n_ahead), drop = FALSE] * scale + center)
  }
  attr(result, "tsp") <- c(min(prd_time), max(prd_time), tsp_train[3])

  result <- list(fcst = result, y = train)
  class(result) <- c("bvarprd", "list")
  return(result)
}
