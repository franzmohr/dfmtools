#' Add the Log Predictive Density of a Forecast
#'
#' Scores the forecast of a dynamic factor model against the observations its
#' horizon realised.
#'
#' @param object an object of class 'dfmodel', usually the result of a call to
#' \code{\link[bvartools]{add_posterior_forecasts}}.
#' @param test_sample a time-series object of the observed series that covers
#' the forecast periods. If \code{NULL} (default), the values in
#' \code{data$test$x} of the object are used, which is what a model carries
#' after it has been scored once.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The log predictive density of period \eqn{T + i} is stored in
#' \code{posterior$forecast$loglik}, one row per draw and one column per scored
#' period. Each column conditions on the observations realised before it, so the
#' log of the mean of the draws of a column is the one step ahead predictive
#' density given everything known up to that period, and those sum over the
#' periods to the log predictive likelihood of the whole realised stretch.
#'
#' A factor model is scored by filtering rather than by evaluating a likelihood
#' on another sample. Its history reaches the density through the factors, which
#' are latent, and its in-sample log-likelihood conditions on the factors the
#' sampler drew -- of which there are none outside the sample. So at every
#' scored period the realised observation updates the distribution of the
#' factors before the next period is predicted, and the column is that period's
#' prediction error decomposition. The filter starts from the drawn factors at
#' the end of the sample and is certain of them: conditional on a draw they are
#' part of what was drawn, and averaging over the draws integrates them out with
#' everything else the posterior carries.
#'
#' Where the loadings, the transition or the volatilities move with time, they
#' take a step of their random walks per scored period, in the order
#' \code{\link[bvartools]{add_posterior_forecasts}} steps them and subject to the same
#' \code{model$forecast_states}. Only the free elements of the loading matrix
#' walk: the identifying block is fixed, and a score that let it move would
#' change the rotation and the scale the factors were estimated under.
#'
#' Fewer realised periods than the horizon is not an error. The periods that are
#' there are the ones that can be scored, and the rest of the forecast is left
#' alone.
#'
#' @return The object in \code{object} with \code{posterior$forecast$loglik}
#' added, a \code{\link[coda]{mcmc}} object with one row per draw and one column
#' per scored period, and with the values it was scored against in
#' \code{data$test$x}.
#'
#' @examples
#'
#' set.seed(1234)
#' data <- matrix(stats::rnorm(200), 50, 4)
#' data <- stats::ts(data, start = c(2000, 1), frequency = 4)
#'
#' train <- stats::window(data, end = c(2011, 1))
#' model <- create_dfmodel(train, n = 1, p = 1, iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model)
#' model <- add_posterior_coefficients(add_initial_values(model))
#' model <- add_posterior_forecasts(model, n_ahead = 3)
#'
#' model <- add_predictive_loglik(model, test_sample = stats::window(data, start = c(2011, 2)))
#' dim(model[["posterior"]][["forecast"]][["loglik"]])
#'
#' @export
#' @method add_predictive_loglik dfmodel
add_predictive_loglik.dfmodel <- function(object, test_sample = NULL, ...) {

  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use 'add_posterior_coefficients' first.")
  }
  if (is.null(object[["posterior"]][["factors"]][["coeffs"]])) {
    stop("Argument 'object' does not contain posterior draws of the factors, which a dynamic ",
         "factor model is scored by filtering on from.")
  }
  if (is.null(object[["model"]][["h"]]) || object[["model"]][["h"]] < 1) {
    stop("Argument 'object' has no forecast horizon. Use 'add_posterior_forecasts' first.")
  }

  realised <- if (is.null(test_sample)) {
    .realised_series(object)
  } else {
    .align_test_series(object, test_sample)
  }
  if (is.null(realised)) {
    return(object)
  }

  class_of_object <- class(object)

  # What a model was scored against travels with it, so that the same model can
  # be scored again without the sample being supplied a second time.
  object[["data"]][["test"]][["x"]] <- realised

  algorithm <- object[["model"]][["algorithm"]]
  if (is.null(algorithm)) {
    stop("Element 'model$algorithm' is missing. Was the object produced by create_dfmodel?")
  }

  if (algorithm == "DfmNormalGamma") {
    object <- .DfmNormalGammaScore(object)
  } else if (algorithm == "DfmNormalStochvol") {
    object <- .DfmNormalStochvolScore(object)
  } else if (algorithm == "DfmTvpGamma") {
    object <- .DfmTvpGammaScore(object)
  } else if (algorithm == "DfmTvpStochvol") {
    object <- .DfmTvpStochvolScore(object)
  } else {
    stop("Algorithm '", algorithm, "' cannot be scored.")
  }

  object[["posterior"]][["forecast"]][["loglik"]] <-
    .mcmc_draws(object[["model"]], object[["posterior"]][["forecast"]][["loglik"]])

  class(object) <- class_of_object

  return(object)
}


# The periods of a test sample that a model's forecast covers, one row per
# period and one column per series, or NULL where the sample does not reach
# them -- which is not an error: a model estimated to the end of a series
# forecasts past what was ever observed.
#
# Matched by time where both the sample and the test data carry one, and taken
# from the top otherwise: a matrix of realised periods is the other way a
# caller has of saying which periods these are.
.align_test_series <- function(object, test_sample) {

  m <- object[["model"]][["m"]]
  h <- object[["model"]][["h"]]
  tsp_train <- stats::tsp(object[["data"]][["x"]])
  tsp_test <- stats::tsp(test_sample)

  test_sample <- stats::na.omit(as.matrix(test_sample))
  if (NCOL(test_sample) != m) {
    stop("Argument 'test_sample' has ", NCOL(test_sample), " columns, but the model has ", m,
         " observed series.")
  }

  if (!is.null(tsp_train) && !is.null(tsp_test)) {
    starts_at <- tsp_train[2] + 1 / tsp_train[3]
    periods <- seq(from = tsp_test[1], by = 1 / tsp_test[3], length.out = nrow(test_sample))
    keep <- periods >= starts_at - 1e-8
    if (!any(keep)) {
      return(NULL)
    }
    test_sample <- test_sample[keep, , drop = FALSE]
  }

  if (nrow(test_sample) == 0) {
    return(NULL)
  }

  return(test_sample[seq_len(min(h, nrow(test_sample))), , drop = FALSE])
}


# The realised values a model carries in data$test$x, checked against what it
# forecast. Already aligned: they are the periods of the horizon and nothing
# else.
.realised_series <- function(object) {

  x <- object[["data"]][["test"]][["x"]]
  if (is.null(x)) {
    stop("No test sample was given and the object carries none in data$test$x. Pass ",
         "'test_sample', or score a model that has been scored once before.")
  }

  x <- as.matrix(x)
  m <- object[["model"]][["m"]]
  h <- object[["model"]][["h"]]
  if (ncol(x) != m) {
    stop("Element data$test$x has ", ncol(x), " columns, but the model has ", m,
         " observed series. It is one row per period and one column per series.")
  }
  if (nrow(x) > h) {
    stop("Element data$test$x holds ", nrow(x), " periods, more than the ", h,
         " this model forecasts.")
  }

  return(x)
}
