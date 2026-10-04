#' Expected Size of a Factor Model
#'
#' Calculates how much memory a dynamic factor model or a factor augmented VAR
#' will take once its posterior draws are complete, block by block, from its
#' specification alone, so that it can be called before a long estimation is
#' started.
#'
#' @param object an object of class 'dfmodel' or 'favarmodel', usually a result
#' of a call to \code{\link{create_dfmodel}} or \code{\link{create_favarmodel}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details Every block of posterior draws is a matrix of numbers with one row
#' per kept draw, and each number takes eight bytes. How many draws are kept is
#' set by \code{iterations}, how many columns each block has by the
#' specification: the number of series \eqn{M}, of factors \eqn{N} and of lags
#' \eqn{p}. \strong{The factors are drawn for every period, and time-varying
#' parameters and stochastic volatility multiply the loadings, the transition
#' or the error precisions by the number of periods as well}, which is how a
#' model of a large panel comes to need gigabytes.
#'
#' Forecast draws are counted once \code{\link{add_forecast_input.dfmodel}} has set a horizon,
#' and the coefficients of deterministic terms where the model has them. The size is what the object will hold in
#' memory; see \code{\link[bvartools]{expected_model_size}} for what that means
#' for the disk space of a file and for the warning
#' \code{\link[bvartools]{add_posterior_coefficients}} gives before it starts.
#'
#' @return A data frame of class 'modelsize', as described on the page of
#' \code{\link[bvartools]{expected_model_size.bvarmodel}}: one row per block,
#' with the columns \code{element}, \code{step}, \code{draws}, \code{columns}
#' and \code{bytes}, and a first row for the data, priors and initial values
#' already in the object. \code{sum(result$bytes)} is the total.
#'
#' @examples
#'
#' data("bem_dfmdata")
#' model <- create_dfmodel(bem_dfmdata, p = 1, n = 1, tvp = TRUE,
#'                         iterations = 5000, burnin = 1000)
#' expected_model_size(model)
#'
#' @export
expected_model_size.dfmodel <- function(object, ...) {

  model <- object[["model"]]
  m <- model[["m"]]
  n <- model[["n"]]
  p <- model[["p"]]
  tt <- NROW(object[["data"]][["x"]])
  draws <- as.numeric(model[["iterations"]])
  sv <- identical(model[["error"]], "sv")
  tvp <- isTRUE(model[["tvp"]])
  periods <- if (tvp) tt else 1
  coef <- "add_posterior_coefficients"
  # Deterministic terms, whose coefficients do not drift.
  n_det <- length(model[["deterministic"]])

  # The loadings are the whole M x N matrix, identifying block included, and a
  # random walk only for the free ones beneath it.
  blocks <- list(
    list("posterior$lambda$coeffs", coef, m * n * periods),
    if (tvp) list("posterior$lambda$sigma", coef, (2 * m - n - 1) * n / 2),
    list("posterior$factors$coeffs", coef, n * tt),
    list("posterior$a$coeffs", coef, n^2 * p * periods),
    if (tvp) list("posterior$a$sigma", coef, n^2 * p),
    list("posterior$u_sigma_inv$coeffs", coef, m * (if (sv) tt else 1)),
    if (sv) list("posterior$u_sigma_inv$sigma", coef, m),
    list("posterior$v_sigma_inv$coeffs", coef, n * (if (sv) tt else 1)),
    if (sv) list("posterior$v_sigma_inv$sigma", coef, n),
    if (n_det > 0) list("posterior$c$coeffs", coef, m * n_det),
    list("posterior$loglik", "add_posterior_loglik", tt),
    if (!is.null(model[["h"]])) list("posterior$forecast$forecasts", "add_posterior_forecasts",
                                     m * model[["h"]]))

  .model_size_table(object, blocks, draws)
}

# The data frame expected_model_size() returns, from a list of blocks, each
# its element, the step that adds it and its number of columns.
.model_size_table <- function(object, blocks, draws) {
  blocks <- Filter(Negate(is.null), blocks)
  columns <- vapply(blocks, function(b) as.numeric(b[[3]]), numeric(1))
  result <- data.frame(
    element = c("data, priors and initial values", vapply(blocks, `[[`, character(1), 1)),
    step = c("", vapply(blocks, `[[`, character(1), 2)),
    draws = c(NA, rep(draws, length(blocks))),
    columns = c(NA, columns),
    bytes = c(as.numeric(utils::object.size(object[names(object) != "posterior"])),
              8 * draws * columns),
    stringsAsFactors = FALSE)
  class(result) <- c("modelsize", "data.frame")
  result
}
