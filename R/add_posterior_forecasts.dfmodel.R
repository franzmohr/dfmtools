#' Forecasting Dynamic Factor Models
#'
#' Produces posterior draws of forecasts of the observable variables of a
#' dynamic factor model.
#'
#' @param object an object of class 'dfmodel', usually, a result of a call to
#' \code{\link{add_posterior_coefficients.dfmodel}} and \code{\link{add_forecast_input.dfmodel}}.
#' @param n_ahead deprecated. The forecast horizon is set by
#' \code{\link{add_forecast_input.dfmodel}}, as for a bvartools model; given here, it is passed on
#' to that function with a warning.
#' @param forecast_states character, what a model whose loadings, transition or
#' volatilities drift does with them over the forecast horizon. \code{"simulate"}
#' carries each draw's random walks forward, one step per period -- the free loadings
#' only, the identifying block staying fixed -- so that the forecasts are draws from
#' the posterior predictive distribution of the estimated model. \code{"hold"} keeps
#' them at their values in the last sample period, which gives forecasts conditional
#' on no further drift and narrower intervals, and is what earlier versions of the
#' package did. If \code{NULL} (default), the value in \code{object$model$forecast_states}
#' is used, and \code{"simulate"} when there is none. A model with constant loadings,
#' transition and volatilities is unaffected.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The forecast takes the route it takes for a bvartools model: the horizon and the
#' out-of-sample data are set by \code{\link{add_forecast_input.dfmodel}}, this function simulates
#' the forecasts, and \code{\link[=predict.dfmodel]{predict}} collects them.
#'
#' A dynamic factor model has no lagged variables among its regressors, so the only out-of-sample
#' data it reads are its deterministic terms, where it has any. Each draw's path is the transition
#' equation run forward from the last \eqn{p} drawn factors,
#' \deqn{f_{T+h} = \sum_{j=1}^{p} A_j f_{T+h-j} + v_{T+h},}
#' with an innovation drawn at every step, and the observable variables read off the
#' loadings, \eqn{x_{T+h} = \lambda f_{T+h} + C d_{T+h} + u_{T+h}}. Both error terms are drawn, so
#' the result is a draw from the posterior predictive distribution rather than a
#' conditional mean.
#'
#' The factor path is therefore required and is what \code{\link{add_posterior_coefficients.dfmodel}}
#' stores in element \code{factors}; the forecast cannot be obtained from the
#' parameters alone.
#'
#' Simulating the states forward steps each random walk by the variance of its
#' innovations: \code{sigma} of \code{lambda} and \code{a} where the coefficients drift,
#' and of \code{u_sigma_inv} and \code{v_sigma_inv} under stochastic volatility. A
#' stochastic volatility model fitted with an earlier version of the package lacks the
#' latter two and stops with an error unless \code{forecast_states = "hold"}.
#'
#' Earlier versions set the horizon here, with \code{n_ahead}, and forecast ten periods ahead
#' when it was not given. Both still work, through \code{\link{add_forecast_input.dfmodel}} and
#' with a warning.
#'
#' @return An object of class 'dfmodel', with element \code{forecast} added to its
#' \code{posterior}. It holds one row per draw and \eqn{h \times M} columns, the
#' horizons stacked within a row in the variable order of the sample, on the scale the model
#' was estimated on. \code{\link[=predict.dfmodel]{predict}} returns them in the units of the data.
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1,
#'                         iterations = 20, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#'
#' # Add priors and initial values
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01),
#'                     u = list(shape = 5, rate = 4),
#'                     a = list(vinv = .01),
#'                     v = list(shape = 5, rate = 4))
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#'
#' # Add the forecast horizon, then simulate the forecasts
#' model <- add_forecast_input(model, n_ahead = 4)
#' model <- add_posterior_forecasts(model)
#'
#' # Collect them
#' pred <- predict(model)
#'
#' @export
add_posterior_forecasts.dfmodel <- function(object, n_ahead = NULL, forecast_states = NULL, ...){

  if (!is.null(forecast_states)) {
    object[["model"]][["forecast_states"]] <- match.arg(forecast_states, c("simulate", "hold"))
  }

  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use add_posterior_coefficients first.")
  }
  if (is.null(object[["posterior"]][["factors"]][["coeffs"]])) {
    stop("Argument 'object' does not contain posterior draws of the factors, which a dynamic factor model forecasts from.")
  }

  class_of_object <- class(object)

  # The horizon and the deterministic terms of the forecast periods, set by
  # add_forecast_input() -- or here, the old way, with a warning.
  object <- .legacy_forecast_horizon(object, n_ahead)
  .check_forecast_input(object)

  algorithm <- object[["model"]][["algorithm"]]
  if (is.null(algorithm)) {
    stop("Element 'model$algorithm' is missing. Was the object produced by create_dfmodel?")
  }

  if (algorithm == "DfmNormalGamma") {
    object <- .DfmNormalGammaForecasts(object)
  } else if (algorithm == "DfmNormalStochvol") {
    object <- .DfmNormalStochvolForecasts(object)
  } else if (algorithm == "DfmTvpGamma") {
    object <- .DfmTvpGammaForecasts(object)
  } else if (algorithm == "DfmTvpStochvol") {
    object <- .DfmTvpStochvolForecasts(object)
  } else {
    stop("Algorithm '", algorithm, "' not supported.")
  }

  object[["posterior"]][["forecast"]][["forecasts"]] <-
    .mcmc_draws(object[["model"]], object[["posterior"]][["forecast"]][["forecasts"]])

  class(object) <- class_of_object

  return(object)
}
