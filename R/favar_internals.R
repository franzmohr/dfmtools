#' Resolve a variable to its position
#'
#' Shared by the loading restrictions of \code{\link{add_priors.favarmodel}} and
#' the lookups of \code{\link{irf.favarmodel}}, so that a name means the same
#' thing and a wrong one reads the same way in both.
#'
#' @param value the name or index provided by the user.
#' @param names the names it is looked up among.
#' @param argument the argument's name, for the error message.
#'
#' @return An integer position.
#'
#' @noRd
.favar_position <- function(value, names, argument) {

  if (is.null(value)) {
    stop("Argument '", argument, "' must be provided. Available are: ",
         paste(names, collapse = ", "), ".")
  }
  if (length(value) != 1) {
    stop("Argument '", argument, "' must be a single variable.")
  }

  if (is.character(value)) {
    position <- which(names == value)
    if (length(position) == 0) {
      stop("Variable '", value, "' is not available as '", argument,
           "'. Available are: ", paste(names, collapse = ", "), ".")
    }
    # A panel column and an observed one may carry the same name, and which of
    # the two was meant is not something to guess at.
    if (length(position) > 1) {
      stop("Name '", value, "' is not unique among the variables available as '",
           argument, "'. Use its index, ", paste(position, collapse = " or "),
           ", instead.")
    }
    return(position)
  }

  position <- suppressWarnings(as.integer(value))
  if (is.na(position) || position < 1 || position > length(names)) {
    stop("Argument '", argument, "' must be between 1 and ", length(names), ".")
  }

  return(position)
}

#' The panel's column names
#'
#' @param object an object of class \code{'favarmodel'}.
#' @param m the number of panel series.
#'
#' @return A character vector of length \code{m}.
#'
#' @noRd
.favar_panel_names <- function(object, m) {

  names <- colnames(object[["data"]][["x"]])
  if (is.null(names)) {
    names <- paste0("x", seq_len(m))
  }

  return(names)
}

#' Propagate an impact vector through the transition
#'
#' The response of the state to a structural shock is linear in the shock's
#' impact vector, so the path is propagated directly rather than through the
#' matrices \eqn{\Psi_h} that would generate it. Shared by
#' \code{\link{irf.favarmodel}} and \code{\link{fevd.favarmodel}}, which differ
#' only in the impact vector they start from and in what they do with the
#' result.
#'
#' @param impact the \eqn{N_s} impact vector.
#' @param a the draw of the transition coefficients, or \code{NULL} when
#' \code{p = 0}.
#' @param ns the width of the state.
#' @param p the lag order.
#' @param n_ahead the horizon.
#'
#' @return An \eqn{N_s \times (n\_ahead + 1)} matrix, one column per horizon.
#'
#' @noRd
.favar_impact_path <- function(impact, a, ns, p, n_ahead) {

  path <- matrix(0, ns, n_ahead + 1)
  path[, 1] <- impact

  if (p > 0 && n_ahead > 0) {
    for (h in seq_len(n_ahead)) {
      for (j in seq_len(min(h, p))) {
        path[, h + 1] <- path[, h + 1] +
          a[, ((j - 1) * ns + 1):(j * ns), drop = FALSE] %*% path[, h - j + 1]
      }
    }
  }

  return(path)
}

#' The covariance of the state innovations of one draw
#'
#' @param v_sigma_inv the matrix of draws of the precision.
#' @param i the draw.
#' @param ns the width of the state.
#'
#' @return An \eqn{N_s \times N_s} covariance matrix.
#'
#' @noRd
.favar_q <- function(v_sigma_inv, i, ns) {

  q <- solve(matrix(v_sigma_inv[i, ], ns, ns))

  # A drawn precision is symmetric, but its inverse is only so up to the
  # solver's error, and chol() is entitled to refuse over it.
  return((q + t(q)) / 2)
}

#' The impact matrix of a set of structural shocks
#'
#' @param q the covariance of the state innovations.
#' @param type \code{"oir"}, \code{"gir"} or \code{"feir"}.
#' @param order the Cholesky ordering, used by \code{"oir"} alone.
#' @param ns the width of the state.
#'
#' @return An \eqn{N_s \times N_s} matrix whose \eqn{j}-th column is the impact
#' of a one-standard-deviation shock to variable \eqn{j}.
#'
#' @noRd
.favar_impact_matrix <- function(q, type, order, ns) {

  if (type == "feir") {
    # A unit shock to one element, with no attempt to orthogonalise: the
    # covariance is ignored rather than decomposed.
    return(diag(1, ns))
  }

  if (type == "gir") {
    # Pesaran and Shin: shock j moves the others by as much as their estimated
    # covariance with it says they do, which needs no ordering.
    return(q %*% diag(1 / sqrt(diag(q)), ns))
  }

  root <- matrix(0, ns, ns)
  root[order, order] <- t(chol(q[order, order, drop = FALSE]))

  return(root)
}

#' Validate a Cholesky ordering
#'
#' @param order the ordering given by the user, or \code{NULL} for the state's
#' own.
#' @param state_names the names of the state, for the error message.
#' @param ns the width of the state.
#'
#' @return An integer permutation of \code{seq_len(ns)}.
#'
#' @noRd
.favar_order <- function(order, state_names, ns) {

  if (is.null(order)) {
    return(seq_len(ns))
  }

  if (is.character(order)) {
    order <- match(order, state_names)
  }
  order <- suppressWarnings(as.integer(order))

  if (anyNA(order) || !setequal(order, seq_len(ns))) {
    stop("Argument 'order' must be a permutation of the ", ns,
         " elements of the state: ", paste(state_names, collapse = ", "), ".")
  }

  return(order)
}
