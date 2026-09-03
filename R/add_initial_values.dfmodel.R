#' Add Initial Values to a Dynamic Factor Model
#'
#' Adds initial values to a dynamic factor model, which was produced by
#' function \code{\link{create_dfmodel}} in combination with \code{\link{add_priors.dfmodel}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_dfmodel}}.
#' @param method a character specifying the method of how initial values are generated.
#' Defaults to \code{"prior"}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' For argument \code{method} the following specifications are possible:
#' \describe{
#'   \item{\code{"prior"}}{Initial values are drawn from the prior. Not possible for uninformative priors.}
#' }
#'
#' Which elements are added depends on argument \code{error} of \code{\link{create_dfmodel}}.
#' For \code{error = "gamma"} the two error precisions are drawn from their gamma priors and stored
#' as the diagonal matrices \code{uinv} and \code{vinv}. For \code{error = "sv"} they are replaced by
#' a log-volatility path each -- \code{u_h}, \eqn{T \times M}, and \code{v_h}, \eqn{T \times N} --
#' together with the state before the first period, \code{u_h_init} and \code{v_h_init}. Each path
#' starts flat, at a draw from the prior on that initial state; a diffuse \code{v_i} therefore gives
#' a starting volatility far from anything the data support, which costs burn-in rather than
#' correctness. The variance of the log-volatility innovations is not among these: it is
#' \code{state_variance} from \code{\link{add_priors.dfmodel}}, which the sampler reads from the
#' prior and redraws every iteration.
#'
#' For a model created with \code{tvp = TRUE} the two coefficient blocks are paths as well. Each
#' gets three elements rather than one: \code{lambda} and \code{a}, matrices with one column per
#' period; \code{lambda_init} and \code{a_init}, the state of the period before the sample, drawn
#' from the prior on it; and \code{lambda_sigma_inv} and \code{a_sigma_inv}, the precision of the
#' state innovations, drawn from the inverse gamma prior added by
#' \code{\link{add_priors.dfmodel}}. Each path starts flat at its own initial state, which costs
#' burn-in rather than correctness -- the sampler redraws the whole path from the data in its first
#' iteration, and what is here only conditions that first draw.
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model data
#' model <- create_dfmodel(x = bem_dfmdata, p = 1:2, n = 1,
#'                          iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01),
#'                     u = list(shape = 5, rate = 4),
#'                     a = list(vinv = .01),
#'                     v = list(shape = 5, rate = 4))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' @export
add_initial_values.dfmodel <- function(object, method = "prior", ...){

  if (method == "prior") {

    tvp <- isTRUE(object$model$tvp)

    # lambda
    if (tvp) {
      object$initial <- .dfm_rw_initial(object, "lambda", object$priors$lambda,
                                        nrow(object$priors$lambda$vinv))
    } else {
      object$initial$lambda <- chol(object$priors$lambda$vinv) %*% stats::rnorm(nrow(object$priors$lambda$vinv))
    }

    error <- object$model$error
    if (is.null(error)) {
      stop("Element 'model$error' is missing. Was the object produced by create_dfmodel?")
    }

    if (error == "gamma") {

      # U
      sigma_shape <- object$priors$u$shape
      sigma_rate <- 1 / object$priors$u$rate
      object$initial$uinv <- diag(1, object$model$m)
      for (i in 1:object$model$m) {
        object$initial$uinv[i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i], rate = sigma_rate[i])
      }
      rm(list = c("sigma_shape", "sigma_rate"))

      # V
      sigma_shape <- object$priors$v$shape
      sigma_rate <- 1 / object$priors$v$rate
      object$initial$vinv <- diag(1, object$model$n)
      for (i in 1:object$model$n) {
        object$initial$vinv[i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i], rate = sigma_rate[i])
      }
      rm(list = c("sigma_shape", "sigma_rate"))

    } else if (error == "sv") {

      # Both log-volatilities, at their own widths: one per observed series for
      # the measurement equation, one per factor for the transition.
      tt <- nrow(object$data$x)

      u_start <- .dfm_sv_initial_state(object$priors$u, object$model$m)
      object$initial$u_h <- matrix(u_start, nrow = tt, ncol = object$model$m, byrow = TRUE)
      object$initial$u_h_init <- matrix(u_start)

      v_start <- .dfm_sv_initial_state(object$priors$v, object$model$n)
      object$initial$v_h <- matrix(v_start, nrow = tt, ncol = object$model$n, byrow = TRUE)
      object$initial$v_h_init <- matrix(v_start)

    } else {
      stop("Error specification '", error, "' not supported.")
    }

    if (object$model$p > 0) {
      # A
      if (tvp) {
        object$initial <- .dfm_rw_initial(object, "a", object$priors$a,
                                          object$model$n^2 * object$model$p)
      } else {
        object$initial$a <- object$priors$a$mu + chol(object$priors$a$vinv) %*% stats::rnorm(object$model$n^2 * object$model$p)
      }
    }
  }

  return(object)
}

# The starting values one random walk coefficient block needs: the state before
# the sample, a path that starts flat at it, and the precision of the state
# innovations. Written under `<name>`, `<name>_init` and `<name>_sigma_inv`,
# which is what the C++ binding reads.
#
# One function for both blocks, so that the loadings and the transition cannot
# drift apart in what they carry. `k` is the width -- freely estimated loadings
# for lambda, transition coefficients for a -- and is passed rather than derived,
# because only one of the two priors has a `mu` to count.
#
# A flat path costs burn-in rather than correctness: the sampler redraws the
# whole path in its first iteration, from data, and what is here only conditions
# that first draw.
.dfm_rw_initial <- function(object, name, prior, k) {

  tt <- nrow(object$data$x)
  initial <- object$initial

  mu <- if (is.null(prior$mu)) rep(0, k) else as.numeric(prior$mu)

  # backsolve(chol(vinv), z) and not chol(vinv) %*% z: with vinv = R'R the first
  # has covariance (R'R)^-1 = vinv^-1, which is the prior, and the second has
  # covariance vinv, which is its inverse. See .dfm_sv_initial_state().
  state <- mu + backsolve(chol(prior$vinv), stats::rnorm(k))

  # The variance of the state innovations, from its own inverse gamma prior. The
  # sampler is handed the precision and flips it back on the way in, which is the
  # convention every time varying model in this family follows.
  variance <- 1 / stats::rgamma(k, shape = as.numeric(prior$shape),
                                rate = as.numeric(prior$rate))

  initial[[name]] <- matrix(state, nrow = k, ncol = tt)
  initial[[paste0(name, "_init")]] <- matrix(state)
  initial[[paste0(name, "_sigma_inv")]] <- diag(1 / variance, k)

  return(initial)
}

# One draw from the prior on the log-volatility before the sample,
# N(mu, v_inv^-1), which is what a stochastic volatility path starts flat at.
#
# `backsolve(chol(v_inv), z)` and not `chol(v_inv) %*% z`: with v_inv = R'R the
# first has covariance (R'R)^-1 = v_inv^-1, which is the prior, and the second
# has covariance v_inv, which is its inverse. The two differ by a factor of v_i
# in the standard deviation, and since this goes through exp() into a variance
# the difference is not cosmetic.
.dfm_sv_initial_state <- function(prior, k) {
  as.numeric(prior$mu) + backsolve(chol(prior$v_inv), stats::rnorm(k))
}
