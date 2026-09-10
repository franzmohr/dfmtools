#' Impulse Response Function of a Factor Augmented VAR
#'
#' Computes the posterior distribution of the response of one variable to a
#' shock to one element of the state.
#'
#' @param x an object of class \code{'favarmodel'} containing posterior draws.
#' The argument is named \code{x} because that is what the generic in
#' \code{bvartools} names it, and not for the panel, which is \code{x$data$x}.
#' @param impulse a character name or integer index of the \emph{state} element
#' the shock hits. See 'Details' for the names.
#' @param response a character name or integer index of the variable whose
#' response is wanted, among the \eqn{k} panel series followed by the
#' \eqn{n_{obs}} observed ones. See 'Details'.
#' @param n_ahead an integer of the response horizon. Defaults to 5.
#' @param ci a numeric between 0 and 1 of the credible interval. Defaults to
#' 0.95.
#' @param shock a numeric multiple of the structural shock. Defaults to 1, which
#' is one standard deviation of the orthogonalised innovation. Negative values
#' give the mirrored response.
#' @param type the identification. \code{"oir"}, the default, orthogonalises the
#' state innovations with the Cholesky factor of \eqn{Q}; \code{"gir"} gives the
#' generalised response of Pesaran and Shin, which needs no ordering;
#' \code{"feir"} leaves the innovations alone. See 'Details'.
#' @param order the Cholesky ordering, as a character or integer permutation of
#' the state. Defaults to the state's own order, factors before observed
#' variables. Only meaningful for \code{type = "oir"}.
#' @param cumulative logical indicating whether the responses should be
#' accumulated over the horizon. Defaults to \code{FALSE}.
#' @param keep_draws logical indicating whether the full posterior of the
#' response should be returned instead of its quantiles. Defaults to
#' \code{FALSE}.
#' @param ... not used.
#'
#' @details Shocks live in the transition equation, so an impulse is always an
#' element of the state \eqn{s_t = (f_t', y_t')'} -- a factor or an observed
#' variable -- and never a series of the panel, whose own error \eqn{e_t} is
#' idiosyncratic by assumption and has no dynamics to propagate. A response, on
#' the other hand, is any \emph{observable}: the \eqn{k} panel series followed by
#' the \eqn{n_{obs}} observed ones, which is the order
#' \code{\link{add_posterior_forecasts.favarmodel}} returns as well. A factor has
#' no entry of its own because it does not need one -- factor \eqn{i} is panel
#' series \eqn{i}, as below. For an observed variable the response is read off
#' the transition,
#' \deqn{\frac{\partial s_{t + h}}{\partial \varepsilon_{t}} = \Psi_h P e_i,
#' \quad \Psi_h = \sum_{j = 1}^{\min(h, p)} \Phi_j \Psi_{h - j}, \quad
#' \Psi_0 = I,}
#' and for a panel series it is carried through the loadings,
#' \eqn{\lambda_{j \cdot} \Psi_h P e_i}, because that is the only way the panel
#' moves at all.
#'
#' \eqn{P} is what argument \code{type} chooses, and it is where the
#' identification lives. Under \code{"oir"} it is the lower Cholesky factor of
#' \eqn{Q} taken in the order given by \code{order}, so a shock is assumed to
#' leave every element ordered before the impulse unmoved on impact and to move
#' all of those after it. That assumption is not tested by anything in the
#' estimated model: \eqn{Q} is unrestricted by construction, which is what a
#' factor augmented VAR is estimated for, and a different ordering gives a
#' different answer from the same draws. Order the state deliberately. Bernanke,
#' Boivin and Eliasz place the slow-moving factors first and the policy rate
#' last, so that policy responds to the factors within the period and they
#' respond to it only with a lag.
#'
#' Under \code{"gir"}, \eqn{P e_i} is replaced by \eqn{Q e_i / \sqrt{q_{ii}}},
#' the generalised response of Pesaran and Shin. There is no ordering to choose:
#' the shock to element \eqn{i} moves the others by as much as their estimated
#' covariance with it says they do, which is a description of the reduced form
#' rather than a structural claim. It answers a different question from
#' \code{"oir"} rather than answering the same one without an assumption, and
#' the two coincide only for the element ordered first. Under \code{"feir"} the
#' innovations are left alone, \eqn{P = I}, which is a reduced-form response and
#' not a structural one at all.
#'
#' Note what \code{"gir"} does to the slow-moving restriction of
#' \code{\link{add_priors.favarmodel}}. That restriction zeroes a slow series'
#' loading on the observed block, so an orthogonalised shock to the policy rate
#' cannot move it on impact. A generalised shock can and generally does, because
#' it moves the factor innovations by their estimated covariance with the policy
#' innovation and the response then reaches the series through its loadings on
#' the factors. The two are not in conflict -- they identify different shocks --
#' but a generalised response is not the Bernanke, Boivin and Eliasz object, and
#' reading it as one gives back the contemporaneous effect the restriction was
#' imposed to rule out.
#'
#' Names are taken from the columns of \code{x} and \code{y} -- which, for a
#' \code{\link[stats]{ts}} built without column names, are R's own
#' \code{"Series 1"} and so on. A factor is named after the panel column that
#' identifies it, the leading \eqn{n \times n} block of \eqn{\lambda_f} being the
#' identity: factor \eqn{i} \emph{is} panel series \eqn{i} up to idiosyncratic
#' noise, and asking for the response of that series returns the factor's, its
#' row of the loadings being a unit vector.
#'
#' Note the scale. With \code{normalize_x = TRUE}, which is
#' \code{\link{create_favarmodel}}'s default, a panel response is in standard
#' deviations of that series rather than in its units, the loading having been
#' estimated against a standardised column. Multiply by
#' \code{attr(x$data$x, "scaled:scale")[response]} to undo it, or estimate with
#' \code{normalize_x = FALSE}.
#'
#' A model with \code{p = 0} has no transition to propagate anything, so its
#' response is the impact alone and zero from horizon one on.
#'
#' @return A time-series object of class \code{'bvarirf'} with one row per
#' horizon, starting at zero, and the lower bound, median and upper bound of the
#' response in its columns. If \code{keep_draws = TRUE}, an \code{'mcmc'} object
#' with one row per draw and one column per horizon instead.
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
#' # The response of the third panel series to a shock to the factor
#' ir <- irf(model, impulse = 1, response = 3, n_ahead = 12)
#'
#' @references
#'
#' Bernanke, B. S., Boivin, J., & Eliasz, P. (2005). Measuring the effects of
#' monetary policy: A factor-augmented vector autoregressive (FAVAR) approach.
#' \emph{The Quarterly Journal of Economics, 120}(1), 387--422.
#'
#' Luetkepohl, H. (2006). \emph{New introduction to multiple time series
#' analysis} (2nd ed.). Berlin: Springer.
#'
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in
#' linear multivariate models. \emph{Economics Letters, 58}(1), 17--29.
#'
#' @seealso \code{\link{fevd.favarmodel}} for the variance decomposition.
#'
#' @export
irf.favarmodel <- function(x, impulse = NULL, response = NULL, n_ahead = 5,
                           ci = 0.95, shock = 1, type = "oir", order = NULL,
                           cumulative = FALSE, keep_draws = FALSE, ...) {

  # 'x' is the generic's name for the model. Everything below calls it what the
  # rest of the package calls it, so that x can go back to meaning the panel.
  object <- x

  if (is.null(object[["posterior"]])) {
    stop("Argument 'x' does not contain posterior draws. Use ",
         "add_posterior_coefficients first.")
  }
  if (n_ahead < 0) {
    stop("Argument 'n_ahead' must be at least 0.")
  }
  if (!is.numeric(shock) || length(shock) != 1) {
    stop("Argument 'shock' must be a single number.")
  }
  if (ci <= 0 || ci >= 1) {
    stop("Argument 'ci' must be between 0 and 1.")
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
  # A factor is the panel series that identifies it, so it borrows that name.
  state_names <- c(panel_names[seq_len(n)], obs_names)

  # Everything that can respond, in the order
  # add_posterior_forecasts.favarmodel returns it: the panel, then the observed
  # block. The factors need no entry of their own, factor i being panel series i.
  response_names <- c(panel_names, obs_names)

  # Impulse and response ----
  impulse <- .favar_position(impulse, state_names, "impulse")
  response <- .favar_position(response, response_names, "response")

  # The observed block is read off the state directly; the panel is carried
  # through the loadings.
  response_in_state <- response > k
  if (response_in_state) {
    response <- n + (response - k)
  }

  # Identification ----
  if (!type %in% c("oir", "gir", "feir")) {
    stop("Argument 'type' must be one of \"oir\", \"gir\" or \"feir\".")
  }
  if (!is.null(order) && type != "oir") {
    stop("Argument 'order' is only meaningful for type \"oir\". A generalised ",
         "response does not order the state, and a forecast error response ",
         "does not orthogonalise it.")
  }
  order <- .favar_order(order, state_names, ns)

  # Draws ----
  v_sigma_inv <- object[["posterior"]][["v_sigma_inv"]][["coeffs"]]
  if (is.null(v_sigma_inv)) {
    stop("Argument 'x' does not contain posterior draws of the state ",
         "covariance, which the identification is taken from.")
  }
  v_sigma_inv <- as.matrix(v_sigma_inv)
  draws <- nrow(v_sigma_inv)

  a <- NULL
  if (p > 0) {
    a <- object[["posterior"]][["a"]][["coeffs"]]
    if (is.null(a)) {
      stop("Argument 'x' does not contain posterior draws of the transition, ",
           "which the response of a model with p > 0 propagates through.")
    }
    a <- as.matrix(a)
  }

  lambda <- NULL
  if (!response_in_state) {
    lambda <- object[["posterior"]][["lambda"]][["coeffs"]]
    if (is.null(lambda)) {
      stop("Argument 'x' does not contain posterior draws of the loadings, ",
           "which the response of a panel series is carried through.")
    }
    lambda <- as.matrix(lambda)
  }

  # Responses ----
  result <- matrix(NA_real_, draws, n_ahead + 1)

  for (i in seq_len(draws)) {

    impact <- .favar_impact_matrix(.favar_q(v_sigma_inv, i, ns), type, order, ns)
    a_i <- if (p > 0) matrix(a[i, ], ns, ns * p) else NULL

    s <- .favar_impact_path(impact[, impulse] * shock, a_i, ns, p, n_ahead)

    if (response_in_state) {
      result[i, ] <- s[response, ]
    } else {
      lambda_i <- matrix(lambda[i, ], k, ns)
      result[i, ] <- lambda_i[response, , drop = FALSE] %*% s
    }
  }

  if (cumulative) {
    result <- t(apply(result, 1, cumsum))
  }

  if (keep_draws) {
    colnames(result) <- 0:n_ahead
    return(coda::as.mcmc(result))
  }

  ci_low <- (1 - ci) / 2
  result <- stats::ts(t(apply(result, 2, stats::quantile,
                              probs = c(ci_low, 0.5, 1 - ci_low))),
                      start = 0, frequency = 1)

  class(result) <- append("bvarirf", class(result))

  return(result)
}
