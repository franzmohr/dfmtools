#' Log-Likelihood of a Dynamic Factor Model
#'
#' Produces the pointwise log-likelihood of a dynamic factor model, one value per
#' posterior draw and observation.
#'
#' @param object an object of class 'dfmodel', usually, a result of a call to
#' \code{\link{add_posterior_coefficients.dfmodel}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The value produced is the measurement density evaluated at the drawn
#' factors,
#' \deqn{p(x_t | f_t, \lambda, U),}
#' and so is the log-likelihood \emph{conditional} on the factor path rather than the
#' marginal one that integrates the factors out. The distinction matters for what an
#' information criterion computed from it means: conditional WAIC assesses a model
#' that treats the states as parameters, which is the quantity a sampler that draws
#' them can report without a filtering pass per draw. The marginal log-likelihood
#' would need a Kalman filter over the sample for every draw and is a different number.
#'
#' The layout is the one \code{waic} and \code{loo} expect: draws along the rows,
#' observations along the columns.
#'
#' @return An object of class 'dfmodel', with element \code{loglik} added to its
#' \code{posterior}.
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
#' # Obtain posterior draws and the log-likelihood
#' model <- add_posterior_coefficients(model)
#' model <- add_posterior_loglik(model)
#'
#' @export
add_posterior_loglik.dfmodel <- function(object, ...){

  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use add_posterior_coefficients first.")
  }
  if (is.null(object[["posterior"]][["factors"]][["coeffs"]])) {
    stop("Argument 'object' does not contain posterior draws of the factors, which the log-likelihood is conditional on.")
  }

  class_of_object <- class(object)

  algorithm <- object[["model"]][["algorithm"]]
  if (is.null(algorithm)) {
    stop("Element 'model$algorithm' is missing. Was the object produced by create_dfmodel?")
  }

  if (algorithm == "DfmNormalGamma") {
    object <- .DfmNormalGammaLogLik(object)
  } else if (algorithm == "DfmNormalStochvol") {
    object <- .DfmNormalStochvolLogLik(object)
  } else if (algorithm == "DfmTvpGamma") {
    object <- .DfmTvpGammaLogLik(object)
  } else if (algorithm == "DfmTvpStochvol") {
    object <- .DfmTvpStochvolLogLik(object)
  } else {
    stop("Algorithm '", algorithm, "' not supported.")
  }

  class(object) <- class_of_object

  return(object)
}
