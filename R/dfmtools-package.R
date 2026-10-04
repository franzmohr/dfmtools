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
#' @importFrom bvartools add_forecast_input
#' @importFrom bvartools prepare_forecast_input
#' @importFrom bvartools add_forecast_errors
#' @importFrom bvartools get_forecast_errors
#' @importFrom stats predict
#' @importFrom bvartools add_posterior_loglik
#' @importFrom bvartools add_predictive_loglik
#' @importFrom bvartools add_seed
#' @importFrom bvartools irf
#' @importFrom bvartools fevd
#' @importFrom bvartools add_sign_zero_restrictions
#' @importFrom bvartools expected_model_size
#' @importFrom bvartools time_variation_test
## usethis namespace: end
NULL

# Unload the DLL when the package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("dfmtools", libpath)
}
