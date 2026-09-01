#' Posterior Simulation of Model Coefficients
#'
#' Forwards model input to posterior simulation functions for dynamic factor models.
#'
#' @param object an object of class 'dfmodel', usually, a result of a call to
#' \code{\link{create_dfmodel}} in combination with \code{\link{add_priors.dfmodel}}
#' and \code{\link{add_initial_values.dfmodel}}.
#' @param posterior_function the function to be applied to the model in argument \code{object}.
#' If \code{NULL} (default), internal functions are used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The function implements the posterior simulation algorithm for Bayesian
#' dynamic factor models described in Chan et al. (2019). The algorithm is implemented
#' in C++ to reduce calculation time.
#'
#' Five blocks are drawn in turn: the path of the unobserved factors, the factor
#' loadings, the precisions of the idiosyncratic errors and of the factor innovations,
#' and the coefficients of the transition equation. The factors are unobserved, so a
#' whole \eqn{T}-period path of them is part of every draw and is returned alongside
#' the parameters in element \code{factors}. It is needed by
#' \code{\link{add_posterior_forecasts.dfmodel}} and
#' \code{\link{add_posterior_loglik.dfmodel}}, neither of which can recompute it.
#'
#' Only the product \eqn{\lambda f_t} is identified, so the leading \eqn{N \times N}
#' block of the loading matrix is fixed unit lower triangular and is not drawn. Element
#' \code{lambda} nevertheless holds the whole \eqn{M \times N} matrix, those fixed ones
#' and zeros included, so that a draw is reshaped rather than unpacked.
#'
#' @return An object of class 'dfmodel', with element \code{posterior} added. It
#' contains the elements \code{lambda}, \code{factors}, \code{a}, \code{u_sigma_inv}
#' and \code{v_sigma_inv}, each a list with element \code{coeffs} holding an object of
#' class \code{\link[coda]{mcmc}} with one row per draw. Note that the two error blocks
#' are precisions, i.e. the inverses of the variances.
#'
#' Within a row, \code{lambda} is the \eqn{M \times N} loading matrix column by column,
#' so a draw is recovered with \code{matrix(draw, m, n)}, and \code{factors} is the
#' path period by period with all \eqn{N} factors of a period together, so
#' \code{t(matrix(draw, n, tt))} recovers a \eqn{T \times N} path. \code{a} holds the
#' transition coefficients as the \eqn{N \times Np} matrix \eqn{[A_1, \dots, A_p]}
#' column by column.
#'
#' @references
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
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
#' # Add priors
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01),
#'                     u = list(shape = 5, rate = 4),
#'                     a = list(vinv = .01),
#'                     v = list(shape = 5, rate = 4))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#'
#' @export
add_posterior_coefficients.dfmodel <- function(object, posterior_function = NULL, ...){

  class_of_object <- class(object)

  # Copy in case the simulation fails
  model <- object
  if ("posterior" %in% names(model)) {
    model[["posterior"]] <- NULL
  }

  if (is.null(posterior_function)) {
    object <- try(
      {
        algorithm <- object[["model"]][["algorithm"]]

        if (is.null(algorithm)) {
          stop("Element 'model$algorithm' is missing. Was the object produced by create_dfmodel?")
        }

        if (algorithm == "DfmNormalGamma") {
          object <- .DfmNormalGammaCoefficients(object)
        } else {
          stop("Algorithm '", algorithm, "' not supported.")
        }

        for (i in c("lambda", "factors", "a", "u_sigma_inv", "v_sigma_inv")) {
          if (!is.null(object[["posterior"]][[i]][["coeffs"]])) {
            object[["posterior"]][[i]][["coeffs"]] <- coda::as.mcmc(object[["posterior"]][[i]][["coeffs"]])
          }
        }

        object
      }
    )
  } else {
    # Apply own function
    object <- try(posterior_function(object))
  }

  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  }

  class(object) <- class_of_object

  return(object)
}
