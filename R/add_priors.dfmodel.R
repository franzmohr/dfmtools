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
#' Every precision here must be larger than 0, and every shape and rate of a gamma that a
#' starting value is drawn from as well. A flat prior is not a specification this package can
#' carry: \code{\link{add_initial_values.dfmodel}} draws each block from its prior, and a
#' precision of 0 has no variance to draw with, as a gamma of shape 0 has no positive value to
#' return. The bound is checked here, where the argument has a name, rather than where the draw
#' breaks.
#'
#' For a model created with \code{tvp = TRUE} both coefficient blocks follow random walks, and both
#' arguments must contain two further elements. \code{vinv} then describes the state of the period
#' before the sample rather than a coefficient that holds throughout: the sampler integrates that
#' state out of the first period, which takes the inverse of its precision. The pair below describes
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
#'   \item{\code{shape}}{a numeric specifying the prior shape parameter, larger than 0. Default is 5.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter, larger than 0. Default is 4.}
#' }
#'
#' For \code{error = "sv"} the log-volatility of every error term follows a random walk, and both
#' arguments must contain all of the following. The names are those \code{bvartools} uses for the
#' stochastic volatility priors of its VAR and VEC models, so that the two families read alike:
#' \describe{
#'   \item{\code{mu}}{a numeric of the prior mean of the initial state of the log-volatilities.}
#'   \item{\code{v_i}}{a numeric of the prior precision of the initial state of the log-volatilities.}
#'   \item{\code{shape}}{a numeric of the prior shape parameter of the variance of the
#'   log-volatility innovations. May be 0, unlike the shapes above: nothing is drawn from this
#'   gamma to start with, the paths starting at \code{state_variance}.}
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
  #
  # Get model specs to obtain total number of coeffs
  m <- object$model$m
  n <- object$model$n
  p <- object$model$p
  error <- object$model$error

  if (is.null(error)) {
    stop("Element 'model$error' is missing. Was the object produced by create_dfmodel?")
  }

  # Total number of freely estimated coefficients in lambda
  n_lambda <- .dfm_n_lambda(m, n)

  # Total # of estimated coefficients in measurement equation
  n_a <- n * n * p

  # Priors for lambda ----
  #
  # A NULL here is a missing argument rather than a way of leaving a prior out:
  # guarding the check with is.null() and then calling diag(NULL, k) only moved
  # the complaint into diag(), where it no longer named the argument. The bound
  # applies where the block has elements; see .check_coefficient_prior().
  .check_coefficient_prior(lambda, "lambda", positive = n_lambda > 0)
  object$priors$lambda <- list(vinv = diag(lambda$vinv, n_lambda))

  # Priors for Phi ----
  #
  # A model with p = 0 has no transition to put a prior on, so `a` is not looked
  # at at all and need not be given.
  if (n_a > 0) {
    .check_coefficient_prior(a, "a")
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
  #
  # `vinv` has to be positive for the sampler to integrate that state out of the
  # first period, and it is checked above for every model rather than for a time
  # varying one alone -- a starting value is drawn from it either way.
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

# The number of freely estimated elements of the M x N loading matrix.
#
# Only the product lambda f_t is identified, so the leading N x N block of lambda
# is fixed unit lower triangular -- ones on the diagonal, zeros above -- which
# pins both the rotation and the scale of the factors. That leaves min(i, N) free
# elements in row i, and summing over the M rows gives N(2M - N - 1)/2.
#
# The count is zero for the one model that has a single observed series and a
# single factor: the whole of lambda is then the identifying block, the loading
# is 1, and what is left is x_t = f_t + u_t with an AR(p) factor -- an
# unobserved-components model, which the sampler estimates like any other. So a
# zero here is a width to carry through, not a specification to reject, and both
# this file and add_initial_values.dfmodel() build empty blocks rather than
# calling chol() on a 0 x 0 matrix.
#
# create_dfmodel() has already rejected N above M, which is the case that would
# make this negative.
.dfm_n_lambda <- function(m, n) (2 * m - n - 1) * n / 2

# The inverse gamma prior on one of the two error precisions. `k` is the width --
# the number of observed series for u, the number of factors for v -- and `name`
# is what a message calls the argument.
.dfm_gamma_prior <- function(spec, k, name) {

  # Named, like the two builders below, rather than counted: a list of length
  # one used to be reported as a length rather than as the field it was short
  # of.
  missing_fields <- setdiff(c("shape", "rate"), names(spec))
  if (length(missing_fields) > 0) {
    stop("Argument '", name, "' is missing the specification",
         if (length(missing_fields) > 1) "s" else "", " ",
         paste0("'", missing_fields, "'", collapse = ", "), ".")
  }
  # Positive rather than non-negative: add_initial_values() draws the precision
  # from this gamma, and rgamma() returns zero for a shape of zero, which is a
  # singular precision rather than a starting value. An improper prior on an
  # error precision is therefore not a specification this package can carry,
  # and saying so here beats saying it one function later.
  .check_prior_number(spec$shape, paste0(name, "$shape"), minimum = 0, strict = TRUE)
  .check_prior_number(spec$rate, paste0(name, "$rate"), minimum = 0, strict = TRUE)

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

  # Positive, for the reason the error precisions are: .dfm_rw_initial() draws
  # the variance of the state innovations from this gamma, and a shape of zero
  # gives a zero there, whose reciprocal is the infinite precision the sampler
  # is handed.
  .check_prior_number(spec$shape, paste0(name, "$shape"), minimum = 0, strict = TRUE)
  .check_prior_number(spec$rate, paste0(name, "$rate"), minimum = 0, strict = TRUE)

  # rep_len() rather than the scalar directly, because `k` is zero for the
  # loadings of a model with a single observed series and matrix(3, 0) is an
  # error rather than an empty matrix. See .dfm_n_lambda() for why that width is
  # a real specification and not a mistake to reject.
  list(shape = matrix(rep_len(spec$shape, k), nrow = k, ncol = 1),
       rate = matrix(rep_len(spec$rate, k), nrow = k, ncol = 1))
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

  # `mu` is a mean and may be anything finite; the rest are bounded. `shape` is
  # allowed to be zero here, unlike in the two builders above, because nothing
  # draws from this gamma to start with: the log-volatilities start at a draw
  # from N(mu, v_i^-1) and the variance of their innovations starts at
  # `state_variance`, so an improper prior on it is a specification this model
  # can carry.
  .check_prior_number(spec$mu, paste0(name, "$mu"))
  .check_prior_number(spec$v_i, paste0(name, "$v_i"), minimum = 0, strict = TRUE)
  .check_prior_number(spec$shape, paste0(name, "$shape"), minimum = 0)
  .check_prior_number(spec$rate, paste0(name, "$rate"), minimum = 0, strict = TRUE)
  .check_prior_number(spec$state_variance, paste0(name, "$state_variance"),
                      minimum = 0, strict = TRUE)
  # Added inside a logarithm, so a zero here is an infinity that would only show
  # up as a broken draw further down.
  .check_prior_number(spec$offset, paste0(name, "$offset"), minimum = 0, strict = TRUE)

  list(mu = matrix(spec$mu, k),
       v_inv = diag(spec$v_i, k),
       shape = matrix(spec$shape, k),
       rate = matrix(spec$rate, k),
       sigma = matrix(spec$state_variance, k),
       offset = matrix(spec$offset, k))
}

# One field of one prior, checked for what it is before it is checked for where
# it lies. Every prior builder in this package goes through this, for both model
# classes, so that a typed or mis-shaped argument is named here rather than
# turning into an NA that surfaces several steps later -- as a character `vinv`
# did, quietly, by way of diag() and a coercion warning.
#
# The type check has to come first and has to be explicit. `spec$shape < 0` on
# the character "a" compares as strings and is FALSE, so a range check on its
# own passes exactly the argument it was written to catch.
#
# `minimum` is the bound and `strict` says whether the bound itself is allowed.
.check_prior_number <- function(value, name, minimum = -Inf, strict = FALSE) {

  if (is.null(value)) {
    stop("Argument '", name, "' is missing.")
  }
  if (!is.numeric(value) || length(value) != 1 || is.na(value) ||
      !is.finite(value)) {
    stop("Argument '", name, "' must be a single finite number.")
  }
  if (strict && value <= minimum) {
    stop("Argument '", name, "' must be larger than ", minimum, ".")
  }
  if (!strict && value < minimum) {
    stop("Argument '", name, "' must be at least ", minimum, ".")
  }

  return(invisible(value))
}

# The normal prior on one of the two coefficient blocks: a single positive
# precision under `vinv`. `name` is what a message calls the argument.
#
# Positive rather than non-negative, which is what the documentation used to
# promise. The starting value of the block is a draw from N(0, vinv^-1), and a
# precision of zero has no variance to draw with: a dynamic factor model's draw
# goes through chol() and stops with "the leading minor of order 1 is not
# positive", and a factor augmented VAR's goes through rnorm(sd = Inf) and
# returns NA, which reaches the sampler and comes back as a complaint about NaN
# from a function the caller has never heard of. Neither is a usable model, so
# the flat prior is refused where it is written rather than where it breaks.
# `positive` says whether the block has anything in it. A model with one series
# and one factor has no freely estimated loading at all, so its `lambda$vinv` is
# never read and a zero there is a number nothing will be drawn from rather than
# a flat prior on something. The type is checked either way.
.check_coefficient_prior <- function(spec, name, positive = TRUE) {

  if (is.null(spec) || !is.list(spec)) {
    stop("Argument '", name, "' must be a named list with element 'vinv'.")
  }
  .check_prior_number(spec$vinv, paste0(name, "$vinv"), minimum = 0,
                      strict = positive)

  return(invisible(spec))
}
