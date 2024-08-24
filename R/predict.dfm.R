#' Predict Method for Objects of Class dfm
#'
#' Forecasting method for class 'dfm' with credible bands.
#'
#' @param object an object of class 'dfm'.
#' @param n_ahead number of steps ahead at which to predict.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param ... additional arguments.
#'
#' @export
predict.dfm <- function(object, n_ahead = 10, ci = 0.95, ...) {

  tt <- nrow(object[["x"]])
  tsp_temp <- stats::tsp(object[["x"]])
  m <- object[["m"]]
  n <- object[["n"]]
  p <- object[["p"]]
  tot <- n * n * p

  draws <- nrow(object[["factor"]])

  result <- array(NA, c(n_ahead, m, draws))
  for (draw in 1:draws) {

    # Obtain forecasts for state equation
    factor_i <- stats::ts(t(matrix(object[["factor"]][draw, ], n)))
    stats::tsp(factor_i) <- tsp_temp
    temp <- bvartools::bvar(y = factor_i,
                            A = t(matrix(object[["a"]][draw, ], 1)),
                            Sigma = t(matrix(object[["v"]][draw, ], 1)))
    states <- predict(temp, n_ahead = n_ahead, ci = ci)[["fcst"]]
    states <- do.call("rbind", lapply(states, function(x) {x[, 1]}))

    # Obtain forecasts for measurement equation
    result[,, draw] <- t(t(matrix(object[["lambda"]][draw, ], n)) %*% states)
    chol_u <- chol(diag(object[["u"]][draw, ], m))
    for (i in 1:n_ahead) {
      result[i,, draw] <- result[i,, draw] + chol_u %*% stats::rnorm(m)
    }
  }

  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  temp <- apply(result, c(1, 2) , stats::quantile, probs = c(ci_low, .5, ci_high))
  result <- c()
  for (i in 1:m) {
    result <- c(result, list(stats::ts(t(temp[,, i]))))
  }
  names(result) <- dimnames(object[["x"]])[[2]]

  if (!is.null(attr(object[["x"]], "ts_info"))) {
    ts_info <- attr(object[["x"]], "ts_info")
    object[["x"]] <- stats::ts(object[["x"]], start = ts_info[1], frequency = ts_info[3])
    attr(object[["x"]], "ts_info") <- NULL

    ts_temp <- stats::ts(0:n_ahead, start = ts_info[2], frequency = ts_info[3])
    ts_temp <- stats::time(ts_temp)[-1]
    for (i in 1:m) {
      stats::tsp(result[[i]]) <- c(ts_temp[1], ts_temp[length(ts_temp)], ts_info[3])
    }
  } else {
    ts_temp <- stats::tsp(object[["x"]])
    for (i in 1:m) {
      result[[i]] <- stats::ts(result[[i]], start = ts_temp[2] + 1 / ts_temp[3], frequency = ts_temp[3])
    }
  }

  result <- list("x" = object[["x"]],
                 "fcst" = result)

  class(result) <- c("dfmprd", "list")
  return(result)
}
