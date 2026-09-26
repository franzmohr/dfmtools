#' @rdname expected_model_size.dfmodel
#' @export
expected_model_size.favarmodel <- function(object, ...) {

  model <- object[["model"]]
  m <- model[["m"]]
  n <- model[["n"]]
  p <- model[["p"]]
  # The VAR runs over the factors and the observed series together.
  k <- n + model[["n_obs"]]
  tt <- NROW(object[["data"]][["x"]])
  draws <- as.numeric(model[["iterations"]])
  coef <- "add_posterior_coefficients"

  blocks <- list(
    list("posterior$lambda$coeffs", coef, m * k),
    list("posterior$factors$coeffs", coef, n * tt),
    list("posterior$a$coeffs", coef, k^2 * p),
    list("posterior$u_sigma_inv$coeffs", coef, m),
    list("posterior$v_sigma_inv$coeffs", coef, k^2),
    list("posterior$loglik", "add_posterior_loglik", tt),
    if (!is.null(model[["h"]])) list("posterior$forecast$forecasts", "add_posterior_forecasts",
                                     (m + model[["n_obs"]]) * model[["h"]]))

  .model_size_table(object, blocks, draws)
}
