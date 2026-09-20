#' Forecast a Factor Augmented VAR
#'
#' Simulates one forecast path per posterior draw.
#'
#' @param object an object of class \code{'favarmodel'} containing posterior
#' draws.
#' @param n_ahead an integer of the forecast horizon. Defaults to 10.
#' @param ... not used.
#'
#' @details There is no forecast design matrix to supply, because a factor
#' augmented VAR has no regressors. Each draw's path is the transition run
#' forward from the last \eqn{p} states,
#' \deqn{s_{T+h} = \sum_{j=1}^{p} \Phi_j s_{T+h-j} + v_{T+h},}
#' with an innovation drawn at every step, and the panel read off the loadings.
#'
#' The horizon covers the \emph{whole} state, observed block included. That is
#' the difference from a dynamic factor model's forecast and it is the point of
#' the model: past the end of the sample the observed variables are no longer
#' data, so they are forecast alongside the factors and their innovations
#' correlate with the factors' through \eqn{Q}.
#'
#' The state path is therefore required, and is what
#' \code{\link{add_posterior_coefficients.favarmodel}} leaves in
#' \code{posterior$factors}.
#'
#' There is no \code{add_predictive_loglik} method for this class, so a
#' forecast made here is not scored against what its horizon realised the way
#' \code{\link{add_predictive_loglik.dfmodel}} scores a dynamic factor model's;
#' see there for why. The draws returned here are what a score would be computed
#' from.
#'
#' @return The model object with element \code{forecast} added to its
#' \code{posterior}. It holds one row per draw and
#' \eqn{h \times (k + n_{obs})} columns -- \strong{wider than the panel}, because
#' each horizon carries the \eqn{k} panel series followed by the
#' \eqn{n_{obs}} observed ones. Those are what a FAVAR is usually forecast for
#' and they have no other object to go in, a dynamic factor model's forecast
#' being \eqn{h \times k}.
#'
#' @examples
#'
#' data("bem_dfmdata")
#'
#' panel <- bem_dfmdata[, -1]
#' observed <- bem_dfmdata[, 1, drop = FALSE]
#'
#' model <- create_favarmodel(x = panel, y = observed, p = 1, n = 1,
#'                            iterations = 500, burnin = 100)
#' # Chosen number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' model <- add_posterior_forecasts(model, n_ahead = 4)
#' dim(model$posterior$forecast$forecasts)
#'
#' @export
add_posterior_forecasts.favarmodel <- function(object, n_ahead = 10, ...) {

  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use ",
         "add_posterior_coefficients first.")
  }
  if (is.null(object[["posterior"]][["factors"]][["coeffs"]])) {
    stop("Argument 'object' does not contain posterior draws of the factors, ",
         "which a factor augmented VAR forecasts from.")
  }
  .check_whole_number(n_ahead, "n_ahead", 1)

  class_of_object <- class(object)

  # The horizon is the whole of what the sampler needs; there is no forecast
  # design matrix for a model without regressors.
  object[["model"]][["h"]] <- as.integer(n_ahead)

  algorithm <- object[["model"]][["algorithm"]]
  if (is.null(algorithm)) {
    stop("Element 'model$algorithm' is missing. Was the object produced by ",
         "create_favarmodel?")
  }
  if (algorithm != "FavarNormalWishart") {
    stop("Algorithm '", algorithm, "' not supported.")
  }

  object <- .FavarNormalWishartForecasts(object)
  object[["posterior"]][["forecast"]][["forecasts"]] <-
    .mcmc_draws(object[["model"]], object[["posterior"]][["forecast"]][["forecasts"]])

  class(object) <- class_of_object

  return(object)
}
