#' Forecast Error Variance Decomposition of a Factor Augmented VAR
#'
#' Decomposes the forecast error variance of one variable into the contributions
#' of the structural shocks to the state.
#'
#' @param x an object of class \code{'favarmodel'} containing posterior draws.
#' The argument is named \code{x} because that is what the generic in
#' \code{bvartools} names it, and not for the panel, which is \code{x$data$x}.
#' @param response a character name or integer index of the variable whose
#' forecast error is decomposed, among the \eqn{k} panel series followed by the
#' \eqn{n_{obs}} observed ones.
#' @param n_ahead an integer of the horizon. Defaults to 5.
#' @param type the identification, \code{"oir"} or \code{"gir"}. See
#' \code{\link{irf.favarmodel}}.
#' @param order the Cholesky ordering, as a character or integer permutation of
#' the state. Only meaningful for \code{type = "oir"}.
#' @param normalise_gir logical indicating whether the shares of a generalised
#' decomposition should be rescaled to account for the whole of the forecast
#' error. Defaults to \code{FALSE}. See 'Details'.
#' @param ... not used.
#'
#' @details The \eqn{h}-step forecast error of the state is
#' \eqn{\sum_{l = 0}^{h} \Psi_l P \varepsilon_{t + h - l}}, so the share of its
#' variance attributable to shock \eqn{j} is
#' \deqn{\omega_{ij}(h) = \frac{\sum_{l = 0}^{h} (e_i' \Psi_l P e_j)^2}
#' {\sum_{l = 0}^{h} e_i' \Psi_l Q \Psi_l' e_i},}
#' with \eqn{\Psi_l} and \eqn{P} as in \code{\link{irf.favarmodel}}. A panel
#' series is reached through its row of the loadings, \eqn{\lambda_{i \cdot}}
#' replacing \eqn{e_i'} throughout.
#'
#' A panel series carries a second source of forecast error, and it is reported
#' rather than hidden. Its idiosyncratic term \eqn{e_{i,t}} is white noise, so
#' it contributes its variance \eqn{\sigma^2_{e,i}} once at every horizon rather
#' than accumulating, but it contributes it to the denominator all the same. The
#' decomposition therefore has one column per state shock and a final column,
#' \code{"idiosyncratic"}, holding the share of the forecast error the common
#' component does not explain. That share is what a factor model claims is
#' specific to the series, and reading it is usually more informative than the
#' shock shares beside it. An observed variable has no such column: it is part of
#' the state and is measured without error, so its forecast error is the state's.
#'
#' Under \code{"oir"} the shares sum to one at every horizon, since
#' \eqn{PP' = Q}. Under \code{"gir"} they do not: the generalised shocks are not
#' orthogonal, so their contributions overlap and sum to more or less than the
#' total. \code{normalise_gir = TRUE} rescales the state shares to fill exactly
#' the share the state does account for, leaving the idiosyncratic column alone,
#' which makes the columns comparable across horizons at the cost of the
#' overlap no longer being visible.
#'
#' Unlike \code{\link{irf.favarmodel}}, which returns the posterior distribution
#' of the response, this returns the posterior \emph{mean} of the shares, one
#' row per horizon -- the convention \code{bvartools} uses, and a reflection of
#' the fact that quantiles of a ratio of variances taken shock by shock need not
#' sum to anything in particular.
#'
#' @return A time-series object of class \code{'bvarfevd'} with one row per
#' horizon, starting at zero, and one column per state element, plus a column
#' \code{"idiosyncratic"} when the response is a panel series.
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' panel <- bem_dfmdata[, -1]
#' observed <- bem_dfmdata[, 1, drop = FALSE]
#'
#' model <- create_favarmodel(x = panel, y = observed, p = 1, n = 1,
#'                            iterations = 500, burnin = 100)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' vd <- fevd(model, response = 3, n_ahead = 12)
#'
#' @references
#'
#' Luetkepohl, H. (2006). \emph{New introduction to multiple time series
#' analysis} (2nd ed.). Berlin: Springer.
#'
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in
#' linear multivariate models. \emph{Economics Letters, 58}(1), 17--29.
#'
#' @seealso \code{\link{irf.favarmodel}} for the responses themselves.
#'
#' @export
fevd.favarmodel <- function(x, response = NULL, n_ahead = 5, type = "oir",
                            order = NULL, normalise_gir = FALSE, ...) {

  object <- x

  if (is.null(object[["posterior"]])) {
    stop("Argument 'x' does not contain posterior draws. Use ",
         "add_posterior_coefficients first.")
  }
  if (n_ahead < 0) {
    stop("Argument 'n_ahead' must be at least 0.")
  }
  if (!type %in% c("oir", "gir")) {
    stop("Argument 'type' must be one of \"oir\" or \"gir\". A decomposition ",
         "of type \"feir\" is not defined: unorthogonalised innovations have ",
         "no separate variances to apportion.")
  }
  if (!is.null(order) && type != "oir") {
    stop("Argument 'order' is only meaningful for type \"oir\".")
  }

  n <- object[["model"]][["n"]]
  n_obs <- object[["model"]][["n_obs"]]
  p <- object[["model"]][["p"]]
  k <- object[["model"]][["m"]]
  ns <- n + n_obs

  # Names ----
  panel_names <- .favar_panel_names(object, k)
  obs_names <- colnames(object[["data"]][["y"]])
  if (is.null(obs_names)) {
    obs_names <- paste0("y", seq_len(n_obs))
  }
  state_names <- c(panel_names[seq_len(n)], obs_names)

  response <- .favar_position(response, c(panel_names, obs_names), "response")
  response_in_state <- response > k
  if (response_in_state) {
    response <- n + (response - k)
  }

  order <- .favar_order(order, state_names, ns)

  # Draws ----
  v_sigma_inv <- object[["posterior"]][["v_sigma_inv"]][["coeffs"]]
  if (is.null(v_sigma_inv)) {
    stop("Argument 'x' does not contain posterior draws of the state ",
         "covariance, which the decomposition is taken from.")
  }
  v_sigma_inv <- as.matrix(v_sigma_inv)
  draws <- nrow(v_sigma_inv)

  a <- NULL
  if (p > 0) {
    a <- object[["posterior"]][["a"]][["coeffs"]]
    if (is.null(a)) {
      stop("Argument 'x' does not contain posterior draws of the transition.")
    }
    a <- as.matrix(a)
  }

  lambda <- u_sigma_inv <- NULL
  if (!response_in_state) {
    lambda <- object[["posterior"]][["lambda"]][["coeffs"]]
    u_sigma_inv <- object[["posterior"]][["u_sigma_inv"]][["coeffs"]]
    if (is.null(lambda) || is.null(u_sigma_inv)) {
      stop("Argument 'x' does not contain posterior draws of the loadings and ",
           "idiosyncratic precisions, which the forecast error of a panel ",
           "series is built from.")
    }
    lambda <- as.matrix(lambda)
    u_sigma_inv <- as.matrix(u_sigma_inv)
  }

  # Shares ----
  width <- ns + if (response_in_state) 0L else 1L
  total <- matrix(0, n_ahead + 1, width)

  for (i in seq_len(draws)) {

    q <- .favar_q(v_sigma_inv, i, ns)
    a_i <- if (p > 0) matrix(a[i, ], ns, ns * p) else NULL

    # The selector that turns a state path into the response: a unit vector for
    # an observed variable, a row of the loadings for a panel series.
    if (response_in_state) {
      selector <- diag(1, ns)[response, ]
      sigma2_e <- 0
    } else {
      selector <- matrix(lambda[i, ], k, ns)[response, ]
      sigma2_e <- 1 / u_sigma_inv[i, response]
    }

    # The denominator is the true forecast error variance and does not depend on
    # which decomposition is reported, so any factor of Q serves for it.
    root <- .favar_impact_matrix(q, "oir", seq_len(ns), ns)
    impact <- if (type == "oir") {
      .favar_impact_matrix(q, "oir", order, ns)
    } else {
      .favar_impact_matrix(q, "gir", order, ns)
    }

    contribution <- numerator <- matrix(0, n_ahead + 1, ns)
    for (j in seq_len(ns)) {
      contribution[, j] <-
        selector %*% .favar_impact_path(root[, j], a_i, ns, p, n_ahead)
      numerator[, j] <-
        selector %*% .favar_impact_path(impact[, j], a_i, ns, p, n_ahead)
    }

    # Accumulated over the horizon, plus the idiosyncratic term, which is white
    # noise and so enters once rather than accumulating.
    state_mse <- cumsum(rowSums(contribution^2))
    denominator <- state_mse + sigma2_e

    shares <- apply(numerator^2, 2, cumsum)
    if (n_ahead == 0) {
      shares <- matrix(shares, 1, ns)
    }
    shares <- shares / denominator

    if (type == "gir" && normalise_gir) {
      overlap <- rowSums(shares)
      keep <- state_mse / denominator
      shares <- shares * ifelse(overlap > 0, keep / overlap, 0)
    }

    if (!response_in_state) {
      shares <- cbind(shares, sigma2_e / denominator)
    }

    total <- total + shares
  }

  result <- total / draws
  colnames(result) <- c(state_names,
                        if (response_in_state) NULL else "idiosyncratic")

  result <- stats::ts(result, start = 0, frequency = 1)
  class(result) <- append("bvarfevd", class(result))

  return(result)
}
