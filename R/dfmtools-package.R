#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib dfmtools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import methods
#' @importFrom bvartools add_priors
#' @importFrom bvartools add_initial_values
#' @importFrom bvartools add_posterior_coefficients
#' @importFrom bvartools add_posterior_forecasts
#' @importFrom bvartools add_posterior_loglik
## usethis namespace: end
NULL

# Unload the DLL when the package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("dfmtools", libpath)
}
