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
#'   \item{\code{omega_v}}{a positive numeric, in place of \code{shape} and \code{rate}: the
#'   variance of the normal prior on the signed standard deviation of the state innovations. See
#'   'The non-centred prior' below.}
#' }
#'
#' @section The non-centred prior:
#'
#' A random walk \eqn{x_t = x_{t-1} + v_t}, \eqn{v_t \sim N(0, \sigma)}, is the same model as
#' \eqn{x_t = x_0 + \omega \tilde{x}_t} with \eqn{\tilde{x}_t} a standard random walk and
#' \eqn{\omega = \pm\sqrt{\sigma}}, and \code{omega_v} is the variance of the normal prior
#' \eqn{\omega \sim N(0, V_\omega)} on that signed standard deviation. It replaces \code{shape}
#' and \code{rate} rather than joining them -- one random walk takes one prior on how far it moves
#' -- and a group that gives both is refused. The implied prior on the variance itself is
#' \eqn{\mathrm{Gamma}(1/2, 1 / (2 V_\omega))}, which puts more mass near zero than the inverse
#' gamma does, so a block that does not move is held at rest more readily.
#'
#' What it buys is a test. \eqn{\omega = 0} is a constant coefficient and an interior point of
#' this prior's support, where a variance of zero is the boundary of the gamma's, so the
#' Savage-Dickey density ratio of Chan (2018) is available:
#' \code{\link[=time_variation_test.dfmodel]{time_variation_test}} reads it off the draws and
#' reports a Bayes factor per state and per block.
#'
#' It is available for the random walks of the two samplers that draw them that way, and for no
#' other block: \code{lambda} and \code{a} of a model created with \code{tvp = TRUE}, and, for one
#' created with \code{tvp = TRUE} and \code{error = "sv"}, the log-volatility groups \code{u} and
#' \code{v} as well. The four are independent -- a model may take it for its loadings and leave
#' the transition on the gamma prior. A model with \code{error = "sv"} and constant coefficients
#' has drifting log-volatilities but no non-centred draw for them, so \code{omega_v} is refused
#' there as everywhere else it would not be read.
#'
#' The stochastic volatility groups keep their remaining elements under this prior:
#' \code{state_variance} is still the starting value of the variance of the log-volatility
#' innovations, and \code{mu}, \code{v_i} and \code{offset} are still read. Only \code{shape} and
#' \code{rate} are replaced.
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
#'   \item{\code{omega_v}}{a positive numeric, in place of \code{shape} and \code{rate}: the
#'   variance of the normal prior on the signed standard deviation of the log-volatility
#'   innovations. See 'The non-centred prior' below.}
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
#' # The same with the non-centred prior on both random walks, which
#' # time_variation_test() can test for time variation afterwards
#' model_nc <- add_priors(model_tvp,
#'                        lambda = list(vinv = .01, omega_v = .1),
#'                        a = list(vinv = .01, omega_v = .1),
#'                        u = list(shape = 5, rate = 4),
#'                        v = list(shape = 5, rate = 4))
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
  } else {
    # A block that does not drift has no random walk to put this prior on, and
    # nothing downstream would read it, so it is refused rather than stored.
    .dfm_refuse_omega_v(lambda, "lambda", "a model with tvp = TRUE")
    .dfm_refuse_omega_v(a, "a", "a model with tvp = TRUE")
  }

  # Error terms ----
  #
  # The two differ only in their width, so both go through the same builder and
  # neither can end up with a field the other has not got.
  if (error == "gamma") {
    object$priors$u <- .dfm_gamma_prior(u, m, "u")
    object$priors$v <- .dfm_gamma_prior(v, n, "v")
  } else if (error == "sv") {
    # Only DfmTvpStochvol draws a log-volatility non-centred. DfmNormalStochvol
    # takes the same prior group but has no such draw, and the vendored core
    # refuses the group -- by the shape it is then missing, which names neither
    # `omega_v` nor the reason -- so the pair is refused here instead.
    noncentred <- isTRUE(object$model$tvp)
    object$priors$u <- .dfm_sv_prior(u, m, "u", noncentred)
    object$priors$v <- .dfm_sv_prior(v, n, "v", noncentred)
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

  # This error term is a constant precision rather than a path, so there is no
  # random walk here for the non-centred prior to be on.
  .dfm_refuse_omega_v(spec, name, "a model with error = \"sv\" and tvp = TRUE")

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

  # The non-centred parameterisation replaces the pair rather than joining it.
  if (!is.null(spec$omega_v)) {
    return(list(omega_v = .dfm_omega_v_prior(spec, k, name)))
  }

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
.dfm_sv_prior <- function(spec, k, name, noncentred_allowed = TRUE) {

  if (!noncentred_allowed) {
    .dfm_refuse_omega_v(spec, name,
                        "a model with error = \"sv\" and tvp = TRUE")
  }

  # `shape` and `rate` are the prior on the variance of the log-volatility
  # innovations, and `omega_v` is the other way of writing it, so a group that
  # gives the second needs neither of the first. Everything else is required
  # either way: the initial state's prior, the starting value of that variance
  # and the offset are all still read.
  noncentred <- !is.null(spec$omega_v)
  required <- c("mu", "v_i", if (!noncentred) c("shape", "rate"),
                "state_variance", "offset")
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
  if (!noncentred) {
    .check_prior_number(spec$shape, paste0(name, "$shape"), minimum = 0)
    .check_prior_number(spec$rate, paste0(name, "$rate"), minimum = 0, strict = TRUE)
  }
  .check_prior_number(spec$state_variance, paste0(name, "$state_variance"),
                      minimum = 0, strict = TRUE)
  # Added inside a logarithm, so a zero here is an infinity that would only show
  # up as a broken draw further down.
  .check_prior_number(spec$offset, paste0(name, "$offset"), minimum = 0, strict = TRUE)

  prior <- list(mu = matrix(spec$mu, k),
                v_inv = diag(spec$v_i, k),
                sigma = matrix(spec$state_variance, k),
                offset = matrix(spec$offset, k))

  if (noncentred) {
    prior$omega_v <- .dfm_omega_v_prior(spec, k, name)
  } else {
    prior$shape <- matrix(spec$shape, k)
    prior$rate <- matrix(spec$rate, k)
  }

  return(prior)
}

# A group given 'omega_v' where nothing drifts. The vendored core reads that
# prior for the random walk of a block, so a model without one would carry a
# prior that is written, stored and never looked at -- which is the one failure
# mode a prior has no way of announcing later. `where` names the models that do
# have the random walk in question.
.dfm_refuse_omega_v <- function(spec, name, where) {

  if (is.list(spec) && !is.null(spec$omega_v)) {
    stop("Argument '", name, "$omega_v' is only available for ", where, ".",
         call. = FALSE)
  }

  return(invisible(NULL))
}

# The non-centred prior on how far one random walk moves, checked and widened to
# the block. `k` is that block's width and `name` what a message calls the
# argument.
#
# The parameterisation is Fruehwirth-Schnatter and Wagner (2010): the path is
# written x_t = x_0 + omega * xtilde_t with xtilde a standard random walk, so
# `omega_v` is the variance of the normal prior on the signed standard deviation
# omega rather than on the variance omega^2 itself. It replaces `shape` and
# `rate` instead of joining them -- one random walk has one prior on how far it
# moves -- and the vendored core refuses a block that carries both, so the pair
# is refused here, where the argument still has a name.
#
# What it buys is the test: omega = 0 is a constant coefficient and an interior
# point of this prior's support, where a variance of zero is the boundary of the
# gamma's, so the Savage-Dickey density ratio of Chan (2018) can be read off the
# draws. time_variation_test() is what reads it.
.dfm_omega_v_prior <- function(spec, k, name) {

  given <- intersect(c("shape", "rate"), names(spec))
  if (length(given) > 0) {
    stop("Argument '", name, "' gives both 'omega_v' and ",
         paste0("'", given, "'", collapse = ", "),
         ". Use 'omega_v' for the non-centred prior or 'shape' and 'rate' for ",
         "the gamma prior, not both.", call. = FALSE)
  }
  .check_prior_number(spec$omega_v, paste0(name, "$omega_v"), minimum = 0,
                      strict = TRUE)

  # rep_len() for the reason .dfm_rw_prior() gives: `k` is zero for the loadings
  # of a model with a single observed series, and matrix(1, 0) is an error.
  matrix(rep_len(spec$omega_v, k), nrow = k, ncol = 1)
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
