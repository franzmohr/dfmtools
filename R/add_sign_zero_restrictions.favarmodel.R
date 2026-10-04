#' Sign and Zero Restrictions for a Factor Augmented VAR
#'
#' Identifies the shocks of the state equation of an object of class
#' \code{'favarmodel'} by the signs of the impulse responses they produce,
#' together with responses restricted to be exactly zero.
#'
#' @param object an object of class \code{'favarmodel'} containing posterior
#' draws.
#' @param restrictions a data frame of the restrictions, with one row per
#' restriction and the columns \code{impulse}, \code{response}, \code{sign} and
#' optionally \code{horizon}, as in
#' \code{\link[bvartools]{add_sign_zero_restrictions}}. \code{impulse} and
#' \code{response} name elements of the state. See 'Details'.
#' @param draws integer. The number of draws to resample. Defaults to
#' \code{NULL}, which uses the effective sample size of the importance sampler.
#' @param max_tries integer. The largest number of rotations drawn for each
#' posterior draw. Defaults to 1, the algorithm of Arias, Rubio-Ramirez and
#' Waggoner (2018). See \code{\link[bvartools]{add_sign_zero_restrictions}} for
#' why more tries keep the sample exact.
#' @param smooth logical. Should the importance weights be Pareto smoothed?
#' Defaults to \code{TRUE}.
#' @param one_sided logical. Should the numerical derivative behind the
#' importance weights be taken on one side only? Defaults to \code{FALSE}.
#' @param ... not used.
#'
#' @details The shocks of a factor augmented VAR live in its transition
#' equation, the VAR of the state \eqn{s_t = (f_t', y_t')'} with innovation
#' covariance \eqn{Q}. This function identifies them the way
#' \code{\link[bvartools]{add_sign_zero_restrictions}} identifies the shocks of
#' a VAR: a rotation of the Cholesky factor of \eqn{Q} per posterior draw,
#' drawn by the algorithm of Arias, Rubio-Ramirez and Waggoner (2018), which
#' \code{\link[bvartools]{arias_rubio_ramirez_waggoner_2018}} runs on the
#' draws of the transition and of \eqn{Q}. The draws are then resampled by
#' their importance weights, and \code{\link{irf.favarmodel}} and
#' \code{\link{fevd.favarmodel}} read the rotations under \code{type = "sign"}.
#'
#' \strong{Restrictions are on elements of the state}: the \eqn{n} factors and
#' the \eqn{n_{obs}} observed variables. Their names are those
#' \code{\link{irf.favarmodel}} uses for an impulse. A factor is named after
#' the panel series that identifies it -- the leading \eqn{n \times n} block
#' of the loadings is the identity, so factor \eqn{i} is panel series
#' \eqn{i} -- and a restriction on it is a restriction on that series. Each
#' shock is named after a state element as well, and the columns of the
#' rotation are built in the order of the state, so a shock carrying many
#' zero restrictions has to be named after an element early in it.
#'
#' The other panel series cannot be restricted. They respond through their
#' loadings, which differ from draw to draw, while the algorithm restricts
#' the responses through matrices that are the same for every draw -- which is
#' what its importance weights are derived for. A sign on such a series would
#' be something else, and is refused rather than approximated.
#'
#' A table without a zero restriction is accepted: the rotations are then
#' drawn uniformly, the weights are equal, and the result is the posterior
#' under sign restrictions alone.
#'
#' @return The object with its posterior draws resampled, the rotations in
#' element \code{q} of its \code{posterior}, one row per resampled draw, and the
#' restrictions with the diagnostics of the importance sample in element
#' \code{sign_zero_restrictions} of its \code{model}: \code{max_tries}, the
#' number of candidate draws, of rotations drawn (\code{tries}) and of draws
#' that satisfied the signs (\code{accepted}), the \code{effective_sample_size},
#' the \code{max_weight_share} of the largest weight and the \code{pareto_k}
#' shape of the weight tail.
#'
#' @examples
#'
#' data("bem_dfmdata")
#' panel <- bem_dfmdata[, -1]
#' observed <- bem_dfmdata[, 1, drop = FALSE]
#'
#' model <- create_favarmodel(x = panel, y = observed, p = 1, n = 1,
#'                            iterations = 300, burnin = 100)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' # The shock named after the factor raises it on impact and leaves the
#' # observed variable unmoved: a slow-moving factor.
#' state <- c(colnames(panel)[1], colnames(observed))
#' restrictions <- data.frame(impulse = state[1], response = state,
#'                            sign = c(1, 0))
#'
#' set.seed(1)
#' model <- add_sign_zero_restrictions(model, restrictions)
#' ir <- irf(model, impulse = 1, response = 3, type = "sign")
#'
#' @references
#'
#' Arias, J. E., Rubio-Ramirez, J. F., Waggoner, D. F. (2018). Inference based
#' on structural vector autoregressions identified with sign and zero
#' restrictions: Theory and applications. \emph{Econometrica, 86}(2), 685-720.
#'
#' @seealso \code{\link{irf.favarmodel}}, \code{\link{fevd.favarmodel}}.
#'
#' @export
add_sign_zero_restrictions.favarmodel <- function(object, restrictions, draws = NULL,
                                                  max_tries = 1, smooth = TRUE,
                                                  one_sided = FALSE, ...) {

  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use ",
         "add_posterior_coefficients first.")
  }

  n <- object[["model"]][["n"]]
  n_obs <- object[["model"]][["n_obs"]]
  p <- object[["model"]][["p"]]
  k <- object[["model"]][["m"]]
  ns <- n + n_obs
  state_names <- .favar_state_names(object, k, n, n_obs)

  # A panel series beyond the n that identify the factors responds through its
  # loadings, which the algorithm cannot restrict; say so by name rather than
  # letting the lookup call it unknown.
  if (is.data.frame(restrictions)) {
    panel_names <- .favar_panel_names(object, k)
    named <- as.character(c(restrictions[["impulse"]], restrictions[["response"]]))
    beyond <- unique(named[named %in% panel_names[-seq_len(n)] & !named %in% state_names])
    if (length(beyond) > 0) {
      stop("Only elements of the state can be restricted: ",
           paste(state_names, collapse = ", "), ". ",
           paste0("'", beyond, "'", collapse = ", "),
           if (length(beyond) > 1) " are panel series" else " is a panel series",
           " that responds through its loadings, which differ from draw to draw. ",
           "See 'Details'.", call. = FALSE)
    }
  }

  v_sigma_inv <- object[["posterior"]][["v_sigma_inv"]][["coeffs"]]
  if (is.null(v_sigma_inv)) {
    stop("Argument 'object' does not contain posterior draws of the state ",
         "covariance, which the identification rotates.")
  }
  v_sigma_inv <- as.matrix(v_sigma_inv)
  store <- nrow(v_sigma_inv)

  a <- NULL
  if (p > 0) {
    a <- object[["posterior"]][["a"]][["coeffs"]]
    if (is.null(a)) {
      stop("Argument 'object' does not contain posterior draws of the transition.")
    }
    a <- as.matrix(a)
  }

  state_draws <- lapply(seq_len(store), function(i) {
    list(A = if (p > 0) matrix(a[i, ], ns, ns * p) else matrix(0, ns, 0),
         Sigma = .favar_q(v_sigma_inv, i, ns))
  })

  identified <- bvartools::arias_rubio_ramirez_waggoner_2018(
    state_draws, restrictions, state_names, p,
    max_tries = max_tries, smooth = smooth, one_sided = one_sided)

  if (is.null(draws)) {
    draws <- identified[["effective_sample_size"]]
  } else {
    .check_whole_number(draws, "draws", 1)
    draws <- as.integer(draws)
  }

  index <- sample.int(store, size = draws, replace = TRUE, prob = identified[["weights"]])

  object[["posterior"]] <- .favar_resample(object[["posterior"]], index, store)
  object[["posterior"]][["q"]] <- list(
    "coeffs" = coda::mcmc(identified[["q"]][index, , drop = FALSE], start = 1, end = draws,
                          thin = 1)
  )

  object[["model"]][["sign_zero_restrictions"]] <- list(
    "restrictions" = identified[["restrictions"]],
    "max_tries" = as.integer(max_tries),
    "candidates" = store,
    "tries" = identified[["tries"]],
    "accepted" = identified[["accepted"]],
    "effective_sample_size" = identified[["effective_sample_size"]],
    "max_weight_share" = identified[["max_weight_share"]],
    "smooth" = smooth,
    "one_sided" = one_sided,
    "pareto_k" = identified[["pareto_k"]]
  )

  return(object)
}

#' The names of the state
#'
#' A factor is named after the panel series that identifies it, the observed
#' variables after their columns. Shared by the identification and by
#' \code{\link{irf.favarmodel}} and \code{\link{fevd.favarmodel}}, so that an
#' impulse means the same element in all three.
#'
#' @noRd
.favar_state_names <- function(object, k, n, n_obs) {
  panel_names <- .favar_panel_names(object, k)
  obs_names <- colnames(object[["data"]][["y"]])
  if (is.null(obs_names)) {
    obs_names <- paste0("y", seq_len(n_obs))
  }
  c(panel_names[seq_len(n)], obs_names)
}

#' Resample every block of a posterior
#'
#' Anything with one row per draw is a block of draws, and taking some blocks
#' and not others would leave row i belonging to different draws in different
#' blocks. The labels start afresh at one: a resample is not a chain.
#'
#' @noRd
.favar_resample <- function(posterior, index, store) {
  take <- function(x) {
    if ((is.matrix(x) || coda::is.mcmc(x)) && NROW(x) == store) {
      # Not as.matrix(), which names the columns of an unnamed chain "var1",
      # "var2", ... -- names no other draw has, and that a model file does not
      # keep.
      draws <- matrix(unclass(x), nrow = NROW(x), dimnames = dimnames(x))
      return(coda::mcmc(draws[index, , drop = FALSE], start = 1,
                        end = length(index), thin = 1))
    }
    x
  }
  for (i in names(posterior)) {
    element <- posterior[[i]]
    if (is.list(element) && !coda::is.mcmc(element)) {
      for (j in names(element)) {
        element[[j]] <- take(element[[j]])
      }
      posterior[[i]] <- element
    } else {
      posterior[[i]] <- take(element)
    }
  }
  posterior
}

#' The impact matrix of an identified draw
#'
#' The Cholesky factor of \eqn{Q} rotated by the draw's accepted rotation, the
#' convention of \code{\link[bvartools]{arias_rubio_ramirez_waggoner_2018}}.
#'
#' @noRd
.favar_sign_impact <- function(q, rotation, ns) {
  t(chol(q)) %*% matrix(rotation, ns, ns)
}
