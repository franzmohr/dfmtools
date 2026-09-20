#' Log Likelihood of a Factor Augmented VAR
#'
#' Computes the pointwise log likelihood, draws by periods.
#'
#' @param object an object of class \code{'favarmodel'}.
#' @param ... not used.
#'
#' @return The model object with the result attached under \code{posterior}, as
#' element \code{loglik}: a \code{\link[coda]{mcmc}} object with one row per
#' draw and one column per observation.
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
#' model <- add_posterior_loglik(model)
#' dim(model$posterior$loglik)
#'
#' @export
add_posterior_loglik.favarmodel <- function(object, ...) {

  # The same two checks a dynamic factor model makes, for the same reason: the
  # density is conditional on the drawn state path, so a model without one has
  # nothing to evaluate.
  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use ",
         "add_posterior_coefficients first.")
  }
  if (is.null(object[["posterior"]][["factors"]][["coeffs"]])) {
    stop("Argument 'object' does not contain posterior draws of the state, ",
         "which the log-likelihood is conditional on.")
  }

  class_of_object <- class(object)

  algorithm <- object[["model"]][["algorithm"]]
  if (is.null(algorithm)) {
    stop("Element 'model$algorithm' is missing. Was the object produced by ",
         "create_favarmodel?")
  }
  if (algorithm != "FavarNormalWishart") {
    stop("Algorithm '", algorithm, "' not supported.")
  }

  object <- .FavarNormalWishartLogLik(object)

  # A chain like any other, so it carries the labels of the draws that were
  # kept -- the convention the dynamic factor models here already set.
  object[["posterior"]][["loglik"]] <-
    .mcmc_draws(object[["model"]], object[["posterior"]][["loglik"]])

  class(object) <- class_of_object

  return(object)
}
