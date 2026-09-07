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
#' @return The model object with element \code{forecast} added to its
#' \code{posterior}. It holds one row per draw and
#' \eqn{h \times (k + n_{obs})} columns -- **wider than the panel**, because
#' each horizon carries the \eqn{k} panel series followed by the
#' \eqn{n_{obs}} observed ones. Those are what a FAVAR is usually forecast for
#' and they have no other object to go in, a dynamic factor model's forecast
#' being \eqn{h \times k}.
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
  if (n_ahead < 1) {
    stop("Argument 'n_ahead' must be at least 1.")
  }

  class_of_object <- class(object)

  # The horizon is the whole of what the sampler needs; there is no forecast
  # design matrix for a model without regressors.
  object[["model"]][["h"]] <- as.integer(n_ahead)

  model <- object

  object <- try({
    algorithm <- object[["model"]][["algorithm"]]
    if (is.null(algorithm)) {
      stop("Element 'model$algorithm' is missing. Was the object produced by ",
           "create_favarmodel?")
    }
    if (algorithm != "FavarNormalWishart") {
      stop("Algorithm '", algorithm, "' not supported.")
    }

    object <- .FavarNormalWishartForecasts(object)
    object[["posterior"]][["forecast"]] <-
      coda::as.mcmc(object[["posterior"]][["forecast"]])

    object
  })

  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  }

  class(object) <- class_of_object

  return(object)
}
