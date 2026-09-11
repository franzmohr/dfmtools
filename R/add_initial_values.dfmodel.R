#' Add Initial Values to a Dynamic Factor Model
#'
#' Adds initial values to a dynamic factor model, which was produced by
#' function \code{\link{create_dfmodel}} in combination with \code{\link{add_priors.dfmodel}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_dfmodel}}.
#' @param method a character specifying the method of how initial values are generated.
#' Defaults to \code{"pca"}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' For argument \code{method} the following specifications are possible:
#' \describe{
#'   \item{\code{"pca"}}{The loadings start at the estimate the leading \eqn{N} principal
#'   components of \code{x} imply, and every other block is drawn from its prior.}
#'   \item{\code{"prior"}}{Every block is drawn from its prior, the loadings included. Not
#'   possible for uninformative priors, and see below.}
#' }
#'
#' The loadings are the one block whose starting value decides which mode the sampler converges
#' to rather than only how long it takes to get there, which is why they do not begin at a draw.
#' Only the product \eqn{\lambda f_t} is identified, and the restriction that pins it fixes the
#' leading \eqn{N \times N} block of \eqn{\lambda} rather than anything about the factors, so a
#' start whose free loadings are negative where the data want them positive describes a coherent
#' model: the mirror image, in which the factor is the negative of the common component and the
#' series whose loading is fixed at one is treated as noise, its idiosyncratic variance absorbing
#' nearly all of its variation. The sampler converges to that mirror and stays there. It fits far
#' worse, and burn-in does not escape it, since turning the factor around would have to pass
#' through configurations no single Gibbs step will take. Under \code{method = "prior"} the draw
#' comes from a prior whose default precision of 0.01 is a standard deviation of ten, so its signs
#' are close to a coin toss and a sizeable share of seeds end up in the mirror.
#'
#' \code{method = "pca"} has no such freedom. The principal components estimate is what the data
#' say, and it is unique up to the sign convention \code{\link[base]{svd}} happens to return,
#' which the rotation onto the identifying restriction cancels. That rotation -- onto a leading
#' \eqn{N \times N} block equal to the identity, which is the unit lower triangular matrix the
#' restriction asks for with its free elements at zero -- is the one step \code{"pca"} can fail
#' at, and it fails where the first \eqn{N} series hardly load on the leading components at all.
#' Such a panel is worth reordering rather than starting from, so the function says so and falls
#' back to the prior draw.
#'
#' What this removes is the cause rather than the possibility. A start on the right side of the
#' likelihood is not a guarantee on a sample too short or too weakly correlated to hold the chain
#' there, and on a four-series panel of sixty quarters the mirror is still reached occasionally.
#' The sign of the estimated loadings is therefore worth a look whatever the starting values were:
#' a factor of real activity whose loadings are negative where the panel is positively correlated
#' is the symptom, and a longer sample, a broader series in the identifying position or simply
#' another seed is the answer.
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
#' # Or with every block drawn from its prior, the loadings included
#' model <- add_initial_values(model, method = "prior")
#'
#' @references
#'
#' Stock, J. H., & Watson, M. W. (2002). Forecasting using principal components from a large
#' number of predictors. \emph{Journal of the American Statistical Association 97}(460),
#' 1167--1179.
#'
#' @export
add_initial_values.dfmodel <- function(object, method = "pca", ...){

  if (!method %in% c("pca", "prior")) {
    stop("Argument 'method' can be 'pca' or 'prior' for dynamic factor models.")
  }

  tvp <- isTRUE(object$model$tvp)

  # lambda, the block whose starting value decides which mode the sampler
  # ends up in rather than only how long it takes to get there.
  n_lambda <- nrow(object$priors$lambda$vinv)
  lambda <- NULL
  if (method == "pca") {
    lambda <- .dfm_pca_initial(object$data$x, object$model$n, n_lambda)
  }
  if (is.null(lambda)) {
    lambda <- .dfm_normal_initial(0, object$priors$lambda$vinv, n_lambda)
  }

  if (tvp) {
    object$initial <- .dfm_rw_initial(object, "lambda", object$priors$lambda,
                                      n_lambda, state = lambda)
  } else {
    object$initial$lambda <- matrix(lambda, nrow = n_lambda, ncol = 1)
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

  return(object)
}

# The freely estimated loadings that the leading `n` principal components of the
# panel imply, in the order the R side stores them -- column by column, each
# column from the diagonal down -- or NULL where that estimate cannot be rotated
# onto the identifying restriction.
#
# `svd(x)$v` holds the component loadings, and the model wants the leading n x n
# block of them unit lower triangular. Post-multiplying by the inverse of that
# block puts the identity there, which is the unit lower triangular matrix with
# its free elements at zero, and rescales the rest to match. The rotation is also
# what makes the result independent of the arbitrary sign svd() returns each
# column with: flipping a column of v flips the same column of the block, and the
# two cancel.
#
# A block that cannot be inverted is a panel whose first n series hold almost
# none of the common variation. There is nothing to start from there, so the
# caller is told and falls back to the prior.
.dfm_pca_initial <- function(x, n, k) {

  if (k == 0) {
    return(numeric(0))
  }

  m <- ncol(x)
  lambda <- svd(unclass(x), nu = 0, nv = n)$v

  block <- lambda[seq_len(n), , drop = FALSE]
  rotated <- try(lambda %*% solve(block), silent = TRUE)

  if (inherits(rotated, "try-error") || !all(is.finite(rotated))) {
    warning("The first ", n, " series of 'x' carry too little of the common ",
            "variation for the loadings to be started from principal ",
            "components, so they were drawn from the prior instead. Ordering a ",
            "series with a large common component first is the fix.",
            call. = FALSE)
    return(NULL)
  }

  # Column j holds the elements below the diagonal, so a column of the leading
  # block holds fewer of them than a column of the rest.
  unlist(lapply(seq_len(n), function(j) rotated[seq_len(m - j) + j, j]))
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
.dfm_rw_initial <- function(object, name, prior, k, state = NULL) {

  tt <- nrow(object$data$x)
  initial <- object$initial

  # `state` is supplied where the block has a starting value of its own -- the
  # loadings under method = "pca" -- and drawn from the prior on the pre-sample
  # state otherwise.
  if (is.null(state)) {
    mu <- if (is.null(prior$mu)) 0 else prior$mu
    state <- .dfm_normal_initial(mu, prior$vinv, k)
  }

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
