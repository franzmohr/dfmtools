#' Add Priors to a Dynamic Factor Model
#'
#' Adds prior specifications to a list of models, which was produced by
#' function \code{\link{create_dfmodel}}.
#'
#' @param object a list, usually, the output of a call to \code{\link{create_dfmodel}}.
#' @param lambda a named list of prior specifications for the factor loadings in the measurement equation.
#' For the default specification the diagonal elements of the inverse prior variance-covariance matrix are set to 0.01.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' @param u a named list of prior specifications for the error variance-covariance matrix of the
#' measurement equation. Which elements are required depends on argument \code{error} of
#' \code{\link{create_dfmodel}}. See 'Details'.
#' @param a a named list of prior specifications for the coefficients of the transition equation.
#' For the default specification the diagonal elements of the inverse prior variance-covariance matrix are set to 0.01.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' @param v a named list of prior specifications for the error variance-covariance matrix of the
#' transition equation. Which elements are required depends on argument \code{error} of
#' \code{\link{create_dfmodel}}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Argument \code{lambda} can only contain the element \code{vinv}, which is a numeric specifying the prior
#' precision of the loading factors of the measurement equation. Default is 0.01.
#'
#' Argument \code{a} can only contain the element \code{vinv}, which is a numeric specifying the prior
#' precision of the coefficients of the transition equation. Default is 0.01.
#'
#' For a model created with \code{tvp = TRUE} both coefficient blocks follow random walks, and both
#' arguments must contain two further elements. \code{vinv} then describes the state of the period
#' before the sample rather than a coefficient that holds throughout, and the pair below describes
#' how far the state may drift from one period to the next. There are no defaults for them, so a
#' call that forgets is told which elements are missing:
#' \describe{
#'   \item{\code{shape}}{a numeric of the prior shape parameter of the variance of the state
#'   innovations.}
#'   \item{\code{rate}}{a numeric of the prior rate parameter of that variance. The smaller it is,
#'   the more tightly the coefficients are held to a constant.}
#' }
#'
#' The two specifications are independent of one another: a model created with \code{tvp = TRUE}
#' and \code{error = "sv"} takes the state equation above for \code{lambda} and \code{a} and the
#' six-element stochastic volatility specification below for \code{u} and \code{v}. All four groups
#' then carry a \code{shape} and a \code{rate}, at four different widths, and each pair belongs to
#' the random walk of the block it sits in.
#'
#' Arguments \code{u} and \code{v} specify the priors of the two error terms -- \eqn{u_t} of the
#' measurement equation and \eqn{v_t} of the transition equation. Both take the same elements, at
#' different widths: \code{u} describes \eqn{M} observed series and \code{v} describes \eqn{N}
#' factors. Which elements are required depends on how the model was created.
#'
#' For \code{error = "gamma"} the function assumes an inverse gamma prior, and both arguments can
#' contain the following elements:
#' \describe{
#'   \item{\code{shape}}{a numeric specifying the prior shape parameter. Default is 5.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter. Default is 4.}
#' }
#'
#' For \code{error = "sv"} the log-volatility of every error term follows a random walk, and both
#' arguments must contain all of the following. The names are those \code{bvartools} uses for the
#' stochastic volatility priors of its VAR and VEC models, so that the two families read alike:
#' \describe{
#'   \item{\code{mu}}{a numeric of the prior mean of the initial state of the log-volatilities.}
#'   \item{\code{v_i}}{a numeric of the prior precision of the initial state of the log-volatilities.}
#'   \item{\code{shape}}{a numeric of the prior shape parameter of the variance of the
#'   log-volatility innovations.}
#'   \item{\code{rate}}{a numeric of the prior rate parameter of that variance.}
#'   \item{\code{state_variance}}{a numeric of the initial draw for the variance of the
#'   log-volatilities. A starting value rather than a prior; it is kept here because that is where
#'   \code{bvartools} keeps it.}
#'   \item{\code{offset}}{a numeric of the constant added before taking the log of the squared
#'   errors, which keeps that logarithm finite when a residual lands on zero.}
#' }
#'
#' @return A list of models.
#'
#' @references
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#'
#' Omori, Y., Chib, S., Shephard, N., & Nakajima, J. (2007). Stochastic volatility with leverage.
#' Fast and efficient likelihood inference. \emph{Journal of Econometrics 140}(2), 425--449.
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model data
#' model <- create_dfmodel(x = bem_dfmdata, p = 1:2, n = 1,
#'                          iterations = 5000, burnin = 1000)
#' # Number of iterations and burn-in should be much higher.
#'
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01),
#'                     u = list(shape = 5, rate = 4),
#'                     a = list(vinv = .01),
#'                     v = list(shape = 5, rate = 4))
#'
#' # The same for a model with stochastic volatility
#' model_sv <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1, error = "sv",
#'                            iterations = 5000, burnin = 1000)
#' model_sv <- add_priors(model_sv,
#'                        lambda = list(vinv = .01),
#'                        a = list(vinv = .01),
#'                        u = list(mu = 0, v_i = .1, shape = 3, rate = .2,
#'                                 state_variance = .05, offset = 1e-4),
#'                        v = list(mu = 0, v_i = .1, shape = 3, rate = .2,
#'                                 state_variance = .05, offset = 1e-4))
#'
#' # And for a model with time varying coefficients, where both coefficient
#' # blocks take a state equation on top of the prior on their initial state
#' model_tvp <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1, tvp = TRUE,
#'                             iterations = 5000, burnin = 1000)
#' model_tvp <- add_priors(model_tvp,
#'                         lambda = list(vinv = .01, shape = 3, rate = .01),
#'                         a = list(vinv = .01, shape = 3, rate = .01),
#'                         u = list(shape = 5, rate = 4),
#'                         v = list(shape = 5, rate = 4))
#'
#' @export
add_priors.dfmodel <- function(object,
                               lambda = list(vinv = 0.01),
                               u = list(shape = 5, rate = 4),
                               a = list(vinv = 0.01),
                               v = list(shape = 5, rate = 4),
                               ...){

  # Checks - Coefficient priors ----
  if (!is.null(lambda)) {
    if (!is.null(lambda$vinv)) {
      if (lambda$vinv < 0) {
        stop("Argument 'lambda$vinv' must be at least 0.")
      }
    } else {
      stop("Argument 'lambda$vinv' is missing.")
    }
  }

  if (!is.null(a)) {
    if (!is.null(a$vinv)) {
      if (a$vinv < 0) {
        stop("Argument 'a$vinv' must be at least 0.")
      }
    } else {
      stop("Argument 'a$vinv' is missing.")
    }
  }

  # Get model specs to obtain total number of coeffs
  m <- object$model$m
  n <- object$model$n
  p <- object$model$p
  error <- object$model$error

  if (is.null(error)) {
    stop("Element 'model$error' is missing. Was the object produced by create_dfmodel?")
  }

  # Total number of freely estimated coefficients in lambda
  n_lambda <- (2 * m - n - 1) * n / 2

  # Total # of estimated coefficients in measurement equation
  n_a <- n * n * p

  # Priors for lambda ----
  object$priors$lambda <- list(vinv = diag(lambda$vinv, n_lambda))

  # Priors for Phi ----
  if (n_a > 0) {
    object$priors$a <- list(mu = matrix(0, n_a),
                            vinv = diag(a$vinv, n_a))
  }

  # State equations ----
  #
  # Where the coefficients drift, `vinv` above stops being a prior on the
  # coefficients themselves and becomes one on the state of the period before the
  # sample; what is added here is the inverse gamma on the variance of the state
  # innovations. Both blocks take the same pair, at their own widths, and go
  # through the same builder so that neither can end up with a field the other
  # has not got. The names are those bvartools uses for the coefficient prior of
  # its time varying VAR and VEC models.
  if (isTRUE(object$model$tvp)) {
    object$priors$lambda <- c(object$priors$lambda,
                              .dfm_rw_prior(lambda, n_lambda, "lambda"))
    if (n_a > 0) {
      object$priors$a <- c(object$priors$a, .dfm_rw_prior(a, n_a, "a"))
    }
  }

  # Error terms ----
  #
  # The two differ only in their width, so both go through the same builder and
  # neither can end up with a field the other has not got.
  if (error == "gamma") {
    object$priors$u <- .dfm_gamma_prior(u, m, "u")
    object$priors$v <- .dfm_gamma_prior(v, n, "v")
  } else if (error == "sv") {
    object$priors$u <- .dfm_sv_prior(u, m, "u")
    object$priors$v <- .dfm_sv_prior(v, n, "v")
  } else {
    stop("Error specification '", error, "' not supported.")
  }

  return(object)
}

# The inverse gamma prior on one of the two error precisions. `k` is the width --
# the number of observed series for u, the number of factors for v -- and `name`
# is what a message calls the argument.
.dfm_gamma_prior <- function(spec, k, name) {

  if (length(spec) < 2) {
    stop("Argument '", name, "' must be at least of length 2.")
  }
  for (field in c("shape", "rate")) {
    if (!field %in% names(spec)) {
      stop("Argument ", name, "$", field, " is missing.")
    }
  }
  if (spec$shape < 0) {
    stop("Argument '", name, "$shape' must be at least 0.")
  }
  if (spec$rate <= 0) {
    stop("Argument '", name, "$rate' must be larger than 0.")
  }

  list(shape = matrix(spec$shape, k),
       rate = matrix(spec$rate, k))
}

# The state equation of one of the two coefficient blocks: the inverse gamma on
# the variance of the random walk innovations. `k` is the width -- the number of
# freely estimated loadings for lambda, the number of transition coefficients for
# a -- and `name` is what a message calls the argument.
#
# Required rather than defaulted, as the stochastic volatility specification
# below is: a caller who asks for time varying coefficients and supplies no state
# equation has not said how far they may drift, and a default would answer that
# question for them without saying so.
.dfm_rw_prior <- function(spec, k, name) {

  required <- c("shape", "rate")
  missing_fields <- setdiff(required, names(spec))
  if (length(missing_fields) > 0) {
    stop("Argument '", name, "' is missing the state equation ",
         "specification", if (length(missing_fields) > 1) "s" else "", " ",
         paste0("'", missing_fields, "'", collapse = ", "),
         ", which a model with time varying coefficients needs.")
  }

  if (spec$shape < 0) {
    stop("Argument '", name, "$shape' must be at least 0.")
  }
  if (spec$rate <= 0) {
    stop("Argument '", name, "$rate' must be larger than 0.")
  }

  list(shape = matrix(spec$shape, k),
       rate = matrix(spec$rate, k))
}

# The stochastic volatility prior on one of the two error terms. Stored under the
# names bvartools uses for the same object, which is also what the C++ binding
# reads: `sigma` is the starting value of the variance of the log-volatility
# innovations rather than a prior, and `v_inv` is a matrix where the caller gave
# a scalar precision.
.dfm_sv_prior <- function(spec, k, name) {

  required <- c("mu", "v_i", "shape", "rate", "state_variance", "offset")
  missing_fields <- setdiff(required, names(spec))
  if (length(missing_fields) > 0) {
    stop("Argument '", name, "' is missing the stochastic volatility ",
         "specification", if (length(missing_fields) > 1) "s" else "", " ",
         paste0("'", missing_fields, "'", collapse = ", "), ".")
  }

  if (spec$v_i <= 0) {
    stop("Argument '", name, "$v_i' must be larger than 0.")
  }
  if (spec$shape < 0) {
    stop("Argument '", name, "$shape' must be at least 0.")
  }
  if (spec$rate <= 0) {
    stop("Argument '", name, "$rate' must be larger than 0.")
  }
  if (spec$state_variance <= 0) {
    stop("Argument '", name, "$state_variance' must be larger than 0.")
  }
  # Added inside a logarithm, so a zero here is an infinity that would only show
  # up as a broken draw further down.
  if (spec$offset <= 0) {
    stop("Argument '", name, "$offset' must be larger than 0.")
  }

  list(mu = matrix(spec$mu, k),
       v_inv = diag(spec$v_i, k),
       shape = matrix(spec$shape, k),
       rate = matrix(spec$rate, k),
       sigma = matrix(spec$state_variance, k),
       offset = matrix(spec$offset, k))
}
