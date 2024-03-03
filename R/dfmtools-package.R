#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib dfmtools, .registration = TRUE
#' @importFrom coda thin
#' @importFrom Rcpp sourceCpp
#' @import methods
#' @importFrom bvartools add_priors
#' @importFrom bvartools draw_posterior
#' @exportPattern "^[[:alpha:]]+"
## usethis namespace: end
NULL

# Unload the DLL when the package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("dfmtools", libpath)
}
