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
    n_lambda <- nrow(object$priors$lambda$vinv)
    if (tvp) {
      object$initial <- .dfm_rw_initial(object, "lambda", object$priors$lambda, n_lambda)
    } else {
      object$initial$lambda <-
        matrix(.dfm_normal_initial(0, object$priors$lambda$vinv, n_lambda),
               nrow = n_lambda, ncol = 1)
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

      u_start <- .dfm_normal_initial(object$priors$u$mu, object$priors$u$v_inv,
                                     object$model$m)
      object$initial$u_h <- matrix(u_start, nrow = tt, ncol = object$model$m, byrow = TRUE)
      object$initial$u_h_init <- matrix(u_start)

      v_start <- .dfm_normal_initial(object$priors$v$mu, object$priors$v$v_inv,
                                     object$model$n)
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
        n_a <- object$model$n^2 * object$model$p
        object$initial$a <-
          matrix(.dfm_normal_initial(object$priors$a$mu, object$priors$a$vinv, n_a),
                 nrow = n_a, ncol = 1)
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

  mu <- if (is.null(prior$mu)) 0 else prior$mu

  state <- .dfm_normal_initial(mu, prior$vinv, k)

  # The variance of the state innovations, from its own inverse gamma prior. The
  # sampler is handed the precision and flips it back on the way in, which is the
  # convention every time varying model in this family follows. Empty at a width
  # of zero, on the same grounds as the state above.
  variance <- if (k == 0) {
    numeric(0)
  } else {
    1 / stats::rgamma(k, shape = as.numeric(prior$shape),
                      rate = as.numeric(prior$rate))
  }

  initial[[name]] <- matrix(state, nrow = k, ncol = tt)
  initial[[paste0(name, "_init")]] <- matrix(state, nrow = k, ncol = 1)
  initial[[paste0(name, "_sigma_inv")]] <- diag(1 / variance, k)

  return(initial)
}

# One draw of length `k` from a normal prior N(mu, vinv^-1), where `vinv` is a
# precision -- which is what every prior in this file is, whatever it is named.
# Every starting value here that is not a variance goes through this: the
# loadings, the transition, the state before the sample of either when they
# drift, and the log-volatility each error term starts flat at.
#
# `backsolve(chol(vinv), z)` and not `chol(vinv) %*% z`. With vinv = R'R the
# first has covariance (R'R)^-1 = vinv^-1, which is the prior. The second has
# covariance R R', which is neither the prior nor its inverse in general; where
# vinv is diagonal, as every default built by add_priors.dfmodel() is, R R'
# collapses to vinv itself, so the second form draws with standard deviation
# sqrt(vinv) where it should draw with 1/sqrt(vinv).
#
# The two therefore do not differ by a constant, and which way they differ
# inverts with the prior: at the default precision of 0.01 the second form gives
# starting values a hundred times tighter than the prior, and at a precision of
# 100 it gives them a hundred times looser -- so it is the caller who states a
# firm prior who is thrown furthest from it. This file had both forms at once
# until they were reconciled, which is why there is now one function rather than
# a rule to remember at four call sites.
#
# A width of zero is a block with nothing to draw rather than a broken one -- at
# one observed series and one factor no loading is freely estimated, see
# .dfm_n_lambda() -- and chol() has no answer for the 0 x 0 precision that
# leaves, so the empty draw is returned directly.
.dfm_normal_initial <- function(mu, vinv, k) {

  if (k == 0) {
    return(numeric(0))
  }

  as.numeric(mu) + backsolve(chol(vinv), stats::rnorm(k))
}
