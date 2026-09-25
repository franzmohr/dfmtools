#' Test for Time Variation in a Dynamic Factor Model
#'
#' Bayes factors for time variation in the loadings, the factor transition and, under
#' stochastic volatility, the two log-volatility blocks of a model estimated with the
#' non-centred prior \code{omega_v}.
#'
#' @param object an object of class \code{'dfmodel'} with posterior draws, i.e. the output of a
#' call to \code{\link[=add_posterior_coefficients.dfmodel]{add_posterior_coefficients}} for a
#' model whose priors were set with \code{omega_v} in
#' \code{\link[=add_priors.dfmodel]{add_priors}}.
#' @param joint a logical. If \code{TRUE}, the default, a block's table also carries the joint
#' test of "no state in this block moves".
#' @param batches an integer of the number of batches the numerical standard error of each
#' Bayes factor is computed from. At least 2, default is 20.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Under the non-centred parameterisation of Frühwirth-Schnatter and Wagner (2010) a random walk
#' \eqn{x_t = x_{t-1} + v_t}, \eqn{v_t \sim N(0, \sigma)}, is written
#' \eqn{x_t = x_0 + \omega \tilde{x}_t} with \eqn{\tilde{x}_t} a standard random walk and
#' \eqn{\omega = \pm \sqrt{\sigma}}. The prior \eqn{\omega \sim N(0, V_\omega)} that
#' \code{omega_v} sets makes \eqn{\omega = 0} -- the state does not move at all -- an interior
#' point of the support, where a variance of zero is the boundary of the gamma's. The
#' Savage-Dickey density ratio of Chan (2018) can then be read off the draws: the Bayes factor of
#' the time-varying model against the one holding the state at \eqn{x_0} is
#' \eqn{p(\omega = 0) / p(\omega = 0 | y)}, whose denominator the sampler writes out one ordinate
#' per draw.
#'
#' A positive \code{log BF} favours time variation. Values above 3 in absolute terms are usually
#' read as strong evidence, in either direction.
#'
#' The joint row of a block tests "every state in it is constant" against "every one moves",
#' which is not "some state moves": a block in which one loading drifts and many stand still can
#' come out against time variation, each constant state costing about a log point. The per-state
#' rows are what says which is which.
#'
#' Only blocks estimated under \code{omega_v} appear. A model that was given \code{shape} and
#' \code{rate} everywhere carries no ordinates, and the function fails rather than returning an
#' empty table.
#'
#' @return A data frame of class \code{'bvartimevar'} with one row per tested state and columns
#' \code{block}, \code{equation}, \code{term}, \code{log_bf} and \code{nse}, printed by the method
#' \code{bvartools} provides for that class.
#'
#' @references
#'
#' Chan, J. C. C. (2018). Specification tests for time-varying parameter models with stochastic
#' volatility. \emph{Econometric Reviews}, 37(8), 807--823.
#' \doi{10.1080/07474938.2016.1167948}
#'
#' Frühwirth-Schnatter, S., & Wagner, H. (2010). Stochastic model specification search for
#' Gaussian and partial non-Gaussian state space models. \emph{Journal of Econometrics}, 154(1),
#' 85--100. \doi{10.1016/j.jeconom.2009.07.003}
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate a model with time varying coefficients
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1, tvp = TRUE,
#'                         iterations = 20, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#'
#' # Both coefficient blocks under the non-centred prior
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01, omega_v = .1),
#'                     a = list(vinv = .01, omega_v = .1),
#'                     u = list(shape = 5, rate = 4),
#'                     v = list(shape = 5, rate = 4))
#'
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' time_variation_test(model)
#'
#' @export
time_variation_test.dfmodel <- function(object, joint = TRUE, batches = 20, ...) {

  if (!is.logical(joint) || length(joint) != 1 || is.na(joint)) {
    stop("Argument 'joint' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(batches) || length(batches) != 1 || is.na(batches) ||
      batches < 2 || batches != round(batches)) {
    stop("Argument 'batches' must be a single integer of at least 2.", call. = FALSE)
  }

  posterior <- object[["posterior"]]
  if (is.null(posterior)) {
    stop("The model has no posterior draws. Use add_posterior_coefficients() first.",
         call. = FALSE)
  }

  m <- object[["model"]][["m"]]
  n <- object[["model"]][["n"]]
  p <- object[["model"]][["p"]]
  series <- colnames(object[["data"]][["x"]])
  if (is.null(series)) {
    series <- paste0("y", seq_len(m))
  }
  factors <- if (n > 0) paste0("factor", seq_len(n)) else character(0)

  # Where the draws are, which prior group holds the omega_v they are tested
  # against, and what the rows of that block are called. The two volatility
  # blocks exist only under error = "sv"; a model without them simply has no
  # ordinates there, which is the same case as a block left centred.
  blocks <- list(
    list(posterior = "lambda", prior = "lambda", title = "loadings",
         labels = .dfm_time_variation_labels_lambda(m, n, series, factors)),
    list(posterior = "a", prior = "a", title = "transition",
         labels = .dfm_time_variation_labels_a(n, p, factors)),
    list(posterior = "u_sigma_inv", prior = "u", title = "idiosyncratic volatilities",
         labels = list(equation = series, term = rep("log-volatility", m))),
    list(posterior = "v_sigma_inv", prior = "v", title = "factor volatilities",
         labels = list(equation = factors, term = rep("log-volatility", n))))

  result <- NULL
  for (b in blocks) {
    draws <- posterior[[b[["posterior"]]]]
    omega_v <- object[["priors"]][[b[["prior"]]]][["omega_v"]]
    if (is.null(draws[["omega_log_zero"]]) || is.null(omega_v)) {
      next
    }

    log_zero <- .dfm_draws_matrix(draws[["omega_log_zero"]])
    omega_v <- rep(as.numeric(omega_v), length.out = ncol(log_zero))
    prior_zero <- stats::dnorm(0, 0, sqrt(omega_v), log = TRUE)

    equation <- b[["labels"]][["equation"]]
    term <- b[["labels"]][["term"]]
    # A label set that does not fit the draws is dropped rather than recycled
    # into a table that names the wrong state.
    if (length(equation) != ncol(log_zero)) {
      equation <- rep(NA_character_, ncol(log_zero))
      term <- paste0(b[["posterior"]], "[", seq_len(ncol(log_zero)), "]")
    }

    rows <- data.frame(
      block = b[["title"]],
      equation = equation,
      term = term,
      log_bf = prior_zero - apply(log_zero, 2, .dfm_log_mean_exp),
      nse = apply(log_zero, 2, .dfm_nse_log_mean_exp, batches = batches),
      stringsAsFactors = FALSE)

    if (joint && !is.null(draws[["omega_log_zero_joint"]])) {
      log_zero_joint <- as.numeric(.dfm_draws_matrix(draws[["omega_log_zero_joint"]]))
      rows <- rbind(rows, data.frame(
        block = b[["title"]],
        equation = NA_character_,
        term = "(joint)",
        log_bf = sum(prior_zero) - .dfm_log_mean_exp(log_zero_joint),
        nse = .dfm_nse_log_mean_exp(log_zero_joint, batches = batches),
        stringsAsFactors = FALSE))
    }

    result <- rbind(result, rows)
  }

  if (is.null(result)) {
    stop("No block of the model was estimated under the non-centred prior. Set ",
         "'omega_v' in add_priors() for the block to be tested -- 'lambda' or 'a' ",
         "for a model with tvp = TRUE, 'u' or 'v' for one with error = \"sv\" -- ",
         "and draw the posterior again.", call. = FALSE)
  }

  rownames(result) <- NULL
  # The class bvartools prints, so that the two families' tables read alike.
  class(result) <- c("bvartimevar", "data.frame")
  return(result)
}


# The series and the factor of every freely estimated loading, in the order the
# draws are in: R's, which is the column-major order of lower.tri() over the
# M x N loading matrix. The leading N x N block is the identifying one and has
# no free element on or above its diagonal, which is what lower.tri() leaves
# out; every row below it is free across its whole width.
.dfm_time_variation_labels_lambda <- function(m, n, series, factors) {

  if (m < 1 || n < 1) {
    return(list(equation = character(0), term = character(0)))
  }
  free <- which(lower.tri(matrix(0, m, n)), arr.ind = TRUE)
  list(equation = series[free[, "row"]], term = factors[free[, "col"]])
}


# The equation and the regressor of every transition coefficient, in the order
# the draws are in: vec([A_1 ... A_p]), so the factor being explained runs
# fastest and the lagged factor and its lag follow.
.dfm_time_variation_labels_a <- function(n, p, factors) {

  if (is.null(p) || n < 1 || p < 1) {
    return(list(equation = character(0), term = character(0)))
  }
  pos <- seq_len(n * n * p) - 1
  column <- pos %/% n
  list(equation = factors[pos %% n + 1],
       term = paste0(factors[column %% n + 1], ".l", column %/% n + 1))
}


# Draws as a matrix, one row per draw, whatever coda class they arrived in.
.dfm_draws_matrix <- function(draws) {

  attr(draws, "mcpar") <- NULL
  class(draws) <- NULL
  if (is.null(dim(draws))) {
    draws <- matrix(draws, ncol = 1)
  }
  return(draws)
}


# log(mean(exp(x))), shifted by the maximum so that the exponentials of the log
# ordinates -- which are large and negative -- do not underflow to zero before
# they are averaged.
.dfm_log_mean_exp <- function(x) {

  maximum <- max(x)
  maximum + log(mean(exp(x - maximum)))
}


# The numerical standard error of the estimate above, by batch means: the draws
# are cut into `batches` consecutive blocks, and the variance of the block means
# says how much of the estimate is sampling noise rather than posterior. On the
# log scale, so it is comparable with the Bayes factor beside it.
#
# NA where there are fewer than two draws per batch, which is a chain too short
# to say anything about its own error rather than an error to stop for.
.dfm_nse_log_mean_exp <- function(x, batches = 20) {

  batches <- min(batches, floor(length(x) / 2))
  if (batches < 2) {
    return(NA_real_)
  }
  e <- exp(x - max(x))
  means <- vapply(split(e, cut(seq_along(e), batches, labels = FALSE)), mean, numeric(1))
  sqrt(stats::var(means) / batches) / mean(e)
}
