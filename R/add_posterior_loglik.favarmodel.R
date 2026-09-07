#' Log Likelihood of a Factor Augmented VAR
#'
#' Computes the pointwise log likelihood, draws by periods.
#'
#' @param object an object of class \code{'favarmodel'}.
#' @param ... not used.
#'
#' @return The model object with the result attached under \code{posterior}.
#'
#' @export
add_posterior_loglik.favarmodel <- function(object, ...) {

  class_of_object <- class(object)

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

    object <- .FavarNormalWishartLogLik(object)

    object
  })

  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  }

  class(object) <- class_of_object

  return(object)
}
