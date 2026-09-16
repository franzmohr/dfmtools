#' Estimate a Factor Augmented VAR
#'
#' Runs the Gibbs sampler and attaches the posterior draws.
#'
#' @param object an object of class \code{'favarmodel'}.
#' @param posterior_function a function that estimates the model. Used to override the built-in sampler.
#' @param ... not used.
#'
#' @details The internal sampler draws with the seed in
#' \code{object$model$seed}, which \code{\link{add_initial_values.favarmodel}}
#' sets and \code{\link{add_seed.favarmodel}} replaces. R's random number
#' generator is set to that seed, with R's default kinds, for the simulation and
#' put back as it was afterwards. A call of \code{set.seed()} between
#' \code{add_initial_values()} and this function therefore does not change the
#' draws. A model without a seed draws from R's generator as it stands. A
#' \code{posterior_function} is called as it is and decides itself what to do
#' with the seed.
#'
#' @return The model object with the result attached under \code{posterior}.
#'
#' @export
add_posterior_coefficients.favarmodel <- function(object, posterior_function = NULL, ...) {

  class_of_object <- class(object)

  model <- object
  if ("posterior" %in% names(model)) {
    model[["posterior"]] <- NULL
  }

  if (!is.null(posterior_function)) {
    object <- try(posterior_function(object))
    if (inherits(object, "try-error")) {
      object <- c(model, list(error = TRUE))
    }
    class(object) <- class_of_object
    return(object)
  }

  object <- try({
    algorithm <- object[["model"]][["algorithm"]]
    if (is.null(algorithm)) {
      stop("Element 'model$algorithm' is missing. Was the object produced by ",
           "create_favarmodel?")
    }
    if (algorithm != "FavarNormalWishart") {
      stop("Algorithm '", algorithm, "' not supported.")
    }

    object <- .with_model_seed(object[["model"]][["seed"]],
                               .FavarNormalWishartCoefficients(object))

    # Every block is a chain like any other, so every one becomes an mcmc
    # object -- the convention the dynamic factor models here already set.
    for (i in c("lambda", "factors", "a", "u_sigma_inv", "v_sigma_inv")) {
      if (!is.null(object[["posterior"]][[i]][["coeffs"]])) {
        object[["posterior"]][[i]][["coeffs"]] <-
          .mcmc_draws(object[["model"]], object[["posterior"]][[i]][["coeffs"]])
      }
    }

    object
  })

  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  }

  class(object) <- class_of_object

  return(object)
}
