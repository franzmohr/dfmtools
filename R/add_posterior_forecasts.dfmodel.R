#' Forecasting Dynamic Factor Models
#'
#' Produces posterior draws of forecasts of the observable variables of a
#' dynamic factor model.
#'
#' @param object an object of class 'dfmodel', usually, a result of a call to
#' \code{\link{add_posterior_coefficients.dfmodel}}.
#' @param n_ahead an integer of the forecast horizon.
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
#' @details Unlike a VAR or a VEC, a dynamic factor model needs no out-of-sample
#' regressors to forecast from: there are none in the model. Each draw's path is the
#' transition equation run forward from the last \eqn{p} drawn factors,
#' \deqn{f_{T+h} = \sum_{j=1}^{p} A_j f_{T+h-j} + v_{T+h},}
#' with an innovation drawn at every step, and the observable variables read off the
#' loadings, \eqn{x_{T+h} = \lambda f_{T+h} + u_{T+h}}. Both error terms are drawn, so
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
#' @return An object of class 'dfmodel', with element \code{forecast} added to its
#' \code{posterior}. It holds one row per draw and \eqn{h \times M} columns, the
#' horizons stacked within a row in the variable order of the sample.
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
#' # Obtain posterior draws and forecasts
#' model <- add_posterior_coefficients(model)
#' model <- add_posterior_forecasts(model, n_ahead = 4)
#'
#' @export
add_posterior_forecasts.dfmodel <- function(object, n_ahead = 10, forecast_states = NULL, ...){

  if (!is.null(forecast_states)) {
    object[["model"]][["forecast_states"]] <- match.arg(forecast_states, c("simulate", "hold"))
  }

  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use add_posterior_coefficients first.")
  }
  if (is.null(object[["posterior"]][["factors"]][["coeffs"]])) {
    stop("Argument 'object' does not contain posterior draws of the factors, which a dynamic factor model forecasts from.")
  }
  if (n_ahead < 1) {
    stop("Argument 'n_ahead' must be at least 1.")
  }

  class_of_object <- class(object)

  # The horizon is the whole of what the sampler needs; there is no forecast
  # design matrix for a model without regressors.
  object[["model"]][["h"]] <- as.integer(n_ahead)

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
