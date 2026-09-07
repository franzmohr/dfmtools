#' Add Initial Values to a Factor Augmented VAR
#'
#' Adds starting values to an object of class \code{'favarmodel'}.
#'
#' @param object an object of class \code{'favarmodel'}, usually the result of a
#' call to \code{\link{add_priors.favarmodel}}.
#' @param method character specifying how the values are obtained. Only
#' \code{"prior"} is available, which draws each block from its own prior.
#' @param ... not used.
#'
#' @details The loadings arrive as a \eqn{(k - n) \times (n + n_{obs})} matrix --
#' one row per panel series that carries a loading, one column per state element
#' -- rather than as a vector. That is the shape a reader wants, and it removes
#' the one question a vector would raise, which is what order the free elements
#' are in.
#'
#' The precision of the state innovations starts at the mean of its Wishart
#' prior, \eqn{\nu S}, rather than at a draw from it. A draw would be a valid
#' start, but the mean is where a chain spends least time getting away from a
#' badly conditioned first sweep, and \eqn{Q} is the one block here whose
#' conditioning the factor path depends on.
#'
#' @return The model object with an additional element \code{initial}.
#'
#' @examples
#'
#' data("bem_dfmdata")
#'
#' model <- create_favarmodel(x = bem_dfmdata[, -1], y = bem_dfmdata[, 1, drop = FALSE],
#'                            p = 1, n = 1, iterations = 5000, burnin = 1000)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#'
#' @export
add_initial_values.favarmodel <- function(object, method = "prior", ...) {

  if (method != "prior") {
    stop("Only method = 'prior' is available for a factor augmented VAR.")
  }
  if (is.null(object[["priors"]])) {
    stop("Element 'priors' is missing. Did you call add_priors?")
  }

  m <- object[["model"]][["m"]]
  n <- object[["model"]][["n"]]
  n_obs <- object[["model"]][["n_obs"]]
  p <- object[["model"]][["p"]]
  n_state <- n + n_obs

  # Loadings, as the (k - n) x n_state block the binding expects.
  n_lambda <- (m - n) * n_state
  lambda_sd <- sqrt(1 / diag(object[["priors"]][["lambda"]][["vinv"]]))
  object[["initial"]][["lambda"]] <-
    matrix(stats::rnorm(n_lambda, sd = lambda_sd), nrow = m - n, ncol = n_state,
           byrow = TRUE)

  # Transition coefficients.
  if (p > 0) {
    n_a <- n_state * n_state * p
    a_sd <- sqrt(1 / diag(object[["priors"]][["a"]][["vinv"]]))
    object[["initial"]][["a"]] <- matrix(stats::rnorm(n_a, sd = a_sd))
  }

  # Idiosyncratic precisions, one draw per panel series from its own gamma.
  shape <- object[["priors"]][["u"]][["shape"]]
  rate <- object[["priors"]][["u"]][["rate"]]
  object[["initial"]][["uinv"]] <- diag(1, m)
  for (i in seq_len(m)) {
    object[["initial"]][["uinv"]][i, i] <- stats::rgamma(1, shape = shape[i], rate = rate[i])
  }

  # The precision of the state innovations, at the mean of its Wishart prior.
  object[["initial"]][["vinv"]] <-
    object[["priors"]][["v"]][["df"]] * object[["priors"]][["v"]][["scale"]]

  return(object)
}
