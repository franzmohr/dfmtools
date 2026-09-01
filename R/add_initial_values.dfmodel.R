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

    # lambda
    object$initial$lambda <- chol(object$priors$lambda$vinv) %*% stats::rnorm(nrow(object$priors$lambda$vinv))

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
      object$initial$a <- object$priors$a$mu + chol(object$priors$a$vinv) %*% stats::rnorm(object$model$n^2 * object$model$p)
    }
  }

  return(object)
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
