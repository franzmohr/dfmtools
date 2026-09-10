#' Add Priors to a Factor Augmented VAR
#'
#' Adds prior specifications to an object of class \code{'favarmodel'}.
#'
#' @param object an object of class \code{'favarmodel'}, usually the result of a
#' call to \code{\link{create_favarmodel}}.
#' @param lambda a named list of prior specifications for the loadings. Element
#' \code{vinv} is the prior precision of every free element.
#' @param a a named list of prior specifications for the transition
#' coefficients, with element \code{vinv} as above.
#' @param u a named list of prior specifications for the idiosyncratic error
#' precisions, with elements \code{shape} and \code{rate}.
#' @param v a named list of prior specifications for the precision of the state
#' innovations, with elements \code{df} and \code{scale}. Unlike every other
#' error block in this package this one is a full matrix; see 'Details'.
#' @param slow the panel series that do not respond to the observed block within
#' the period, as a character or integer vector, or a named list with elements
#' \code{series} and \code{vinv}. Defaults to \code{NULL}, no restriction. See
#' 'Details'.
#' @param ... not used.
#'
#' @details The idiosyncratic precisions get independent gamma priors, one per
#' panel series, exactly as in \code{\link{add_priors.dfmodel}}. The precision of
#' the state innovations gets a \emph{Wishart} prior instead, because \eqn{Q} is
#' unrestricted in a FAVAR -- its off-diagonal block is the correlation between
#' the factor innovations and the shock to the observed variables, which is the
#' quantity the model exists to measure. \code{df} defaults to the width of the
#' state and \code{scale} to the identity of that size.
#'
#' The number of free loadings is \eqn{(k - n)(n + n_{obs})}: the first \eqn{n}
#' panel series identify the factors and carry no free loading at all, and every
#' series after them loads freely on the whole state. See
#' \code{\link{create_favarmodel}} for why the identification is an identity
#' block rather than the unit lower triangle a dynamic factor model uses.
#'
#' Argument \code{slow} names the series that are assumed not to react to the
#' observed block within the period -- output, employment and prices, against
#' financial variables that reprice on the day. Their loadings on the observed
#' columns of \eqn{\lambda} are pinned at zero, which is the restriction
#' Bernanke, Boivin and Eliasz identify a monetary policy shock with: the
#' factors are then slow-moving too, so a recursive ordering with the policy rate
#' last says something the data has not already been asked to say. Listing one of
#' the first \eqn{n} series is allowed and does nothing, the identification
#' having zeroed their observed columns already -- which is also why those
#' \eqn{n} should themselves be slow-moving series.
#'
#' The restriction is a prior, not a hard zero: element \code{vinv} of the list
#' form sets the precision it is held at, \code{1e12} by default, which leaves a
#' loading at zero to eight decimal places. Nothing stops a series from being
#' slow against one observed variable and not another, but \code{slow} does not
#' express that; write into \code{priors$lambda$vinv} directly for it.
#'
#' @return The model object with an additional element \code{priors}.
#'
#' @examples
#'
#' data("bem_dfmdata")
#'
#' model <- create_favarmodel(x = bem_dfmdata[, -1], y = bem_dfmdata[, 1, drop = FALSE],
#'                            p = 1, n = 1, iterations = 5000, burnin = 1000)
#' model <- add_priors(model)
#'
#' # Real quantities are assumed not to react to the observed block within the
#' # quarter, financial variables to be free to.
#' model <- add_priors(model, slow = c("PCECC96", "PCDGx", "PCESVx", "PCNDx"))
#'
#' @references
#'
#' Bernanke, B. S., Boivin, J., & Eliasz, P. (2005). Measuring the effects of
#' monetary policy: A factor-augmented vector autoregressive (FAVAR) approach.
#' \emph{The Quarterly Journal of Economics, 120}(1), 387--422.
#'
#' @export
add_priors.favarmodel <- function(object,
                                  lambda = list(vinv = 0.01),
                                  a = list(vinv = 0.01),
                                  u = list(shape = 5, rate = 4),
                                  v = list(df = NULL, scale = NULL),
                                  slow = NULL,
                                  ...) {

  m <- object[["model"]][["m"]]
  n <- object[["model"]][["n"]]
  n_obs <- object[["model"]][["n_obs"]]
  p <- object[["model"]][["p"]]

  if (is.null(m) || is.null(n) || is.null(n_obs)) {
    stop("Element 'model' is incomplete. Was the object produced by create_favarmodel?")
  }

  n_state <- n + n_obs

  # The first n panel series identify the factors and carry no free loading;
  # every series after them loads freely on the whole state.
  n_lambda <- (m - n) * n_state
  n_a <- n_state * n_state * p

  if (is.null(lambda[["vinv"]]) || lambda[["vinv"]] < 0) {
    stop("Argument 'lambda$vinv' must be at least 0.")
  }
  object[["priors"]][["lambda"]] <- list(mu = matrix(0, n_lambda),
                                         vinv = diag(lambda[["vinv"]], n_lambda))

  # The slow-moving restriction ----
  if (!is.null(slow)) {

    slow_vinv <- 1e12
    if (is.list(slow)) {
      if (!is.null(slow[["vinv"]])) {
        slow_vinv <- slow[["vinv"]]
      }
      slow <- slow[["series"]]
    }
    if (!is.numeric(slow_vinv) || length(slow_vinv) != 1 || slow_vinv < 0) {
      stop("Argument 'slow$vinv' must be a single number of at least 0.")
    }

    if (length(slow) > 0) {

      series <- vapply(slow, .favar_position, integer(1),
                       names = .favar_panel_names(object, m), argument = "slow")

      # The first n series have no free loading to restrict; the prior mean is
      # zero already, so all that is left is to hold them there.
      series <- unique(series[series > n])

      if (length(series) > 0) {
        position <- as.vector(outer((series - n - 1) * n_state,
                                    n + seq_len(n_obs), "+"))
        object[["priors"]][["lambda"]][["vinv"]][cbind(position, position)] <-
          slow_vinv
      }
    }
  }

  if (n_a > 0) {
    if (is.null(a[["vinv"]]) || a[["vinv"]] < 0) {
      stop("Argument 'a$vinv' must be at least 0.")
    }
    object[["priors"]][["a"]] <- list(mu = matrix(0, n_a),
                                      vinv = diag(a[["vinv"]], n_a))
  }

  for (field in c("shape", "rate")) {
    if (!field %in% names(u)) {
      stop("Argument u$", field, " is missing.")
    }
  }
  if (u[["shape"]] < 0) {
    stop("Argument 'u$shape' must be at least 0.")
  }
  if (u[["rate"]] <= 0) {
    stop("Argument 'u$rate' must be larger than 0.")
  }
  object[["priors"]][["u"]] <- list(shape = matrix(u[["shape"]], m),
                                    rate = matrix(u[["rate"]], m))

  # Q is unrestricted, so its prior is a Wishart rather than a set of
  # independent gammas. The defaults are the least informative proper ones: the
  # width of the state for the degrees of freedom, and the identity for the
  # scale.
  df <- if (is.null(v[["df"]])) n_state else v[["df"]]
  scale <- if (is.null(v[["scale"]])) diag(1, n_state) else v[["scale"]]

  if (df < n_state) {
    stop("Argument 'v$df' must be at least the width of the state (", n_state,
         ") for the Wishart prior to be proper.")
  }
  if (!identical(dim(as.matrix(scale)), as.integer(c(n_state, n_state)))) {
    stop("Argument 'v$scale' must be a ", n_state, " x ", n_state, " matrix.")
  }
  object[["priors"]][["v"]] <- list(df = df, scale = as.matrix(scale))

  return(object)
}
