#' Posterior Simulation
#' 
#' Forwards model input to posterior simulation functions.
#' 
#' @param object a list of model specifications, which should be passed on
#' to function \code{FUN}. Usually, the output of a call to \code{\link{gen_dfm}}
#' in combination with \code{\link{add_priors}}.
#' @param posterior_function the function to be applied to argument \code{object}.
#' If \code{NULL} (default), the internal function \code{\link{dfmpost}} is used.
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of the class of the output of the applied posterior
#' simulation function. In case the package's own function is used, this will
#' result in an object of class \code{"dfm"}.
#' 
#' @examples
#' 
#' # Load data 
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model
#' model <- gen_var(e1, p = 1:2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' @export
draw_posterior.dfmodel <- function(object, posterior_function = NULL, mc.cores = NULL, ...){
  
  if (is.null(posterior_function)) {
    object <- try(dfmpost(object))
  } else {
    # Apply own function
    object <- try(posterior_function(object))
  }
  
  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(object, list(coefficients = NULL))
  }
  
  return(object)
}