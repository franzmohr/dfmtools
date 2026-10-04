#' Prepare Forecast Input for a Factor Model
#'
#' Generates the data a dynamic factor model or a factor augmented VAR is
#' forecast from, the counterpart of \code{\link[bvartools]{prepare_forecast_input}}
#' for a VAR.
#'
#' @param object an object of class \code{'dfmodel'} or \code{'favarmodel'}.
#' @param n_ahead number of steps ahead at which to predict.
#' @param deterministic a time-series object with the deterministic terms of the
#' forecast periods, one column per term of the model. If not specified, the terms are
#' continued from the estimation sample, which works for the constant, the trend and the
#' seasonal dummies that \code{\link{create_dfmodel}} and \code{\link{create_favarmodel}}
#' generate. Only read for a model with deterministic terms.
#' @param ... additional arguments.
#'
#' @details A factor model is forecast by running its transition forward from the drawn
#' factors, so unlike a VAR it has no lagged variables to put into the regressors of a forecast
#' period. What it does need out of sample is the value of its deterministic terms in each
#' forecast period, and that is the whole of \code{x} here. A model without deterministic terms
#' needs nothing but the horizon.
#'
#' A series passed in \code{deterministic} has to cover the \code{n_ahead} periods after the
#' estimation sample at the frequency of the model; one that does not is refused with a message
#' naming the periods it has to cover.
#'
#' @return A list with elements \code{h}, the forecast horizon, and \code{x}, the deterministic
#' terms of the forecast periods: \code{h} rows, one per period, by one column per term, or
#' \code{NULL} for a model without deterministic terms.
#'
#' @examples
#'
#' data("bem_dfmdata")
#'
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1,
#'                         deterministic = "const", seasonal = TRUE,
#'                         iterations = 20, burnin = 10)
#'
#' fcst_input <- prepare_forecast_input(model, n_ahead = 4)
#' fcst_input$x
#'
#' @export
#' @rdname prepare_forecast_input.dfmodel
prepare_forecast_input.dfmodel <- function(object, n_ahead = 10, deterministic = NULL, ...) {
  .check_whole_number(n_ahead, "n_ahead", 1)
  list("h" = as.integer(n_ahead),
       "x" = .forecast_deterministic(object, n_ahead, deterministic))
}

#' @export
#' @rdname prepare_forecast_input.dfmodel
prepare_forecast_input.favarmodel <- function(object, n_ahead = 10, deterministic = NULL, ...) {
  prepare_forecast_input.dfmodel(object, n_ahead = n_ahead, deterministic = deterministic, ...)
}


#' Add Forecast Input to a Factor Model
#'
#' Sets the forecast horizon of a dynamic factor model or a factor augmented VAR and adds the
#' data its forecast is produced from, the counterpart of
#' \code{\link[bvartools]{add_forecast_input}} for a VAR.
#'
#' @param object an object of class \code{'dfmodel'} or \code{'favarmodel'}.
#' @param n_ahead number of steps ahead at which to predict.
#' @param deterministic a time-series object with the deterministic terms of the
#' forecast periods. If not specified, they are continued from the estimation sample. See
#' \code{\link{prepare_forecast_input.dfmodel}}.
#' @param ... arguments passed forward to method.
#'
#' @details This is the first of the three steps a forecast takes, as it is in bvartools:
#' \code{add_forecast_input} sets the horizon and the out-of-sample data,
#' \code{\link{add_posterior_forecasts.dfmodel}} simulates the forecasts, and
#' \code{\link[=predict.dfmodel]{predict}} collects them. The forecast errors of
#' \code{\link{add_forecast_errors.dfmodel}} and the predictive likelihood of
#' \code{\link{add_predictive_loglik.dfmodel}} are computed from what the second step produced.
#'
#' @return The object in \code{object} with \code{model$h} set to \code{n_ahead} and
#' \code{data$forecast$x} added, the deterministic terms of the forecast periods, one row per
#' period; \code{NULL} for a model without deterministic terms.
#'
#' @examples
#'
#' data("bem_dfmdata")
#'
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1,
#'                         deterministic = "const",
#'                         iterations = 20, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' model <- add_forecast_input(model, n_ahead = 4)
#' model <- add_posterior_forecasts(model)
#'
#' @export
#' @rdname add_forecast_input.dfmodel
add_forecast_input.dfmodel <- function(object, n_ahead = 10, deterministic = NULL, ...) {
  fcst_input <- prepare_forecast_input(object, n_ahead = n_ahead,
                                       deterministic = deterministic, ...)
  object[["model"]][["h"]] <- fcst_input[["h"]]
  object[["data"]][["forecast"]] <- list("x" = fcst_input[["x"]])
  return(object)
}

#' @export
#' @rdname add_forecast_input.dfmodel
add_forecast_input.favarmodel <- function(object, n_ahead = 10, deterministic = NULL, ...) {
  add_forecast_input.dfmodel(object, n_ahead = n_ahead, deterministic = deterministic, ...)
}


# What add_posterior_forecasts() does with a horizon given the old way.
#
# Before the forecast took bvartools' route, through add_forecast_input(), it
# was asked for with add_posterior_forecasts(object, n_ahead = ...) -- and an
# object with no horizon at all was forecast ten periods ahead. Both still
# work, through add_forecast_input(), so that a released script runs unchanged,
# and both say what to call instead.
.legacy_forecast_horizon <- function(object, n_ahead) {
  if (!is.null(n_ahead)) {
    warning("Argument 'n_ahead' of add_posterior_forecasts() is deprecated. Set the horizon ",
            "with add_forecast_input() before calling add_posterior_forecasts(), as for a ",
            "bvartools model.", call. = FALSE)
    return(add_forecast_input(object, n_ahead = n_ahead))
  }
  if (is.null(object[["model"]][["h"]])) {
    warning("Model specification does not contain forecast horizon 'h', so it was set to 10. ",
            "Use add_forecast_input() before add_posterior_forecasts(), as for a bvartools ",
            "model; the default will be removed.", call. = FALSE)
    return(add_forecast_input(object, n_ahead = 10))
  }
  object
}


# Refuses a forecast of a model with deterministic terms that was not given
# them for every period of the horizon -- in R, where the message can name the
# function that supplies them, rather than in the sampler.
.check_forecast_input <- function(object) {
  n_det <- length(object[["model"]][["deterministic"]])
  if (n_det == 0) {
    return(invisible(TRUE))
  }
  x <- object[["data"]][["forecast"]][["x"]]
  h <- object[["model"]][["h"]]
  if (is.null(x) || NROW(x) < h || NCOL(x) != n_det) {
    stop("Model specification does not contain the deterministic terms of the ", h,
         " forecast periods. Consider using function add_forecast_input().", call. = FALSE)
  }
  invisible(TRUE)
}
