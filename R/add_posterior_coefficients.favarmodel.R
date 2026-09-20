#' Estimate a Factor Augmented VAR
#'
#' Runs the Gibbs sampler and attaches the posterior draws.
#'
#' @param object an object of class \code{'favarmodel'}.
#' @param posterior_function a function that estimates the model. Used to override the built-in sampler.
#' @param verbose logical indicating whether the sampler should report its
#' progress. Defaults to \code{FALSE}; see
#' \code{\link{add_posterior_coefficients.dfmodel}}.
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
#' A list of models, which \code{\link{create_favarmodel}} returns for more
#' than one lag order or number of factors, is simulated by the list method of
#' \pkg{bvartools}. Its argument \code{cores} shares the models out over that
#' many worker processes, which load this package for the methods they need;
#' see \code{\link[bvartools]{add_posterior_coefficients.modellist}}. Each model
#' draws with its own seed, so the draws do not depend on the number of workers,
#' but the workers run with one thread, and the draws of a factor augmented VAR
#' equal those of a session doing the same rather than of one running several.
#' \code{add_posterior_forecasts()} and \code{add_posterior_loglik()} take
#' \code{cores} for a list as well -- but the guarantee above is about
#' \emph{these} draws and does not extend to them. Only this function is
#' seeded. \code{add_posterior_forecasts()} draws from R's generator as it
#' stands, so its draws depend on the state that generator is in, and on a
#' cluster that is the state of whichever worker took the model: the same list
#' forecast on one worker and on four does not give the same paths.
#' \code{set.seed()} immediately before the call is what fixes them.
#'
#'
#' A simulation that fails does not stop: the model comes back as it went in
#' with an element \code{error} set to \code{TRUE}, so that a list of models
#' loses only the one that failed. Estimation is the only step that does this.
#' \code{add_posterior_forecasts()}, \code{add_posterior_loglik()} and
#' \code{add_predictive_loglik()} are cheap and derived, and what goes wrong in
#' them is nearly always the call rather than the chain, so they report it and
#' stop.
#'
#' @return The model object with the result attached under \code{posterior}.
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
#' dim(model$posterior$lambda$coeffs)
#'
#' @export
add_posterior_coefficients.favarmodel <- function(object, posterior_function = NULL,
                                                  verbose = FALSE, ...) {

  object[["model"]][["verbose"]] <- .check_verbose(verbose)

  class_of_object <- class(object)

  # Both the copy kept for the failure path and the object handed to the
  # sampler are stripped of what a previous run left; see the dynamic factor
  # model's method for why the binding needs that now.
  object[["posterior"]] <- NULL
  object[["error"]] <- NULL
  model <- object

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
