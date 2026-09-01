#' Forecasting Dynamic Factor Models
#'
#' Produces posterior draws of forecasts of the observable variables of a
#' dynamic factor model.
#'
#' @param object an object of class 'dfmodel', usually, a result of a call to
#' \code{\link{add_posterior_coefficients.dfmodel}}.
#' @param n_ahead an integer of the forecast horizon.
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
add_posterior_forecasts.dfmodel <- function(object, n_ahead = 10, ...){

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

  object <- .DfmNormalGammaForecasts(object)

  object[["posterior"]][["forecast"]] <- coda::as.mcmc(object[["posterior"]][["forecast"]])

  class(object) <- class_of_object

  return(object)
}
