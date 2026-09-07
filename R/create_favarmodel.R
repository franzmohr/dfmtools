#' Create a Factor Augmented VAR
#'
#' Produces the input for the estimation of a factor augmented VAR (FAVAR).
#'
#' @param x a time-series object of the stationary panel the factors are
#' extracted from.
#' @param y a time-series object of the observed variables that enter the state
#' vector alongside the factors. These are \emph{not} regressors; see 'Details'.
#' @param p an integer vector of the lag order of the transition equation. See
#' 'Details'.
#' @param n an integer vector of the number of unobserved factors.
#' @param normalize_x logical indicating whether each column of \code{x} should
#' be normalized using \code{scale}. Defaults to \code{TRUE}.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 20000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000).
#'
#' @details The function produces the variable matrices of a factor augmented
#' VAR with measurement equation
#' \deqn{x_t = \lambda_f f_t + \lambda_y y_t + e_t,}
#' where \eqn{x_t} is a \eqn{k \times 1} vector of observed panel series,
#' \eqn{f_t} an \eqn{n \times 1} vector of unobserved factors and \eqn{y_t} the
#' observed variables that also enter the state. \eqn{e_t} is a \eqn{k \times 1}
#' error term with \eqn{e_t \sim N(0, R)} and \eqn{R} diagonal.
#'
#' The transition equation is a VAR in the whole state
#' \eqn{s_t = (f_t', y_t')'},
#' \deqn{s_t = \sum_{i = 1}^{p} \Phi_i s_{t - i} + v_t,}
#' with \eqn{v_t \sim N(0, Q)} and \eqn{Q} unrestricted.
#'
#' **The observed block is part of the state, not a set of regressors.** It
#' appears on the left of the transition equation as well as the right, and the
#' model's own dynamics run over it. That is what a FAVAR is estimated to
#' measure: \eqn{Q} is unrestricted precisely so that its cross block -- the
#' correlation between the factor innovations and the shock to the observed
#' variables -- can be read.
#'
#' Identification differs from a dynamic factor model's, and not by accident.
#' The leading \eqn{n \times n} block of \eqn{\lambda_f} is the \emph{identity}
#' and the observed columns of those same rows are zero, so the first \eqn{n}
#' series of the panel are the factors plus idiosyncratic noise exactly. A
#' dynamic factor model can use a unit lower triangle instead because its
#' \eqn{V} is diagonal, and the two restrictions together admit only the
#' identity rotation. A FAVAR has no diagonal \eqn{Q} to offer -- an unrestricted
#' \eqn{Q} is the model -- so a unit lower triangle on its own would leave the
#' loadings free to wander along a ridge. Order \code{x} so that its first
#' \code{n} columns are the series you are willing to define the factors by.
#'
#' Note what \code{normalize_x} does to the loadings. Each column of \code{x} is
#' divided by its own standard deviation, so an estimated loading is on that
#' standardised scale rather than on the data's -- and the identification pins
#' the factor to the \emph{standardised} first series. That is usually what is
#' wanted, since a factor extracted from series of different scales is a factor
#' of the largest of them. It is also why a loading cannot be compared directly
#' with one implied by the raw data unless the normalisation is undone, or
#' \code{normalize_x = FALSE} is used.
#'
#' If integer vectors are provided as arguments \code{p} or \code{n}, the
#' function produces a distinct model for all combinations of those
#' specifications.
#'
#' @return An object of class \code{'favarmodel'}, which contains the following
#' elements:
#' \item{data}{A list of data objects. Element \code{x} is the normalised panel
#' and element \code{y} the observed factors.}
#' \item{model}{A list of model specifications.}
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # The first series stands in for an observed factor here
#' panel <- bem_dfmdata[, -1]
#' observed <- bem_dfmdata[, 1, drop = FALSE]
#'
#' model <- create_favarmodel(x = panel, y = observed, p = 1, n = 1,
#'                            iterations = 5000, burnin = 1000)
#'
#' @references
#'
#' Bernanke, B. S., Boivin, J., & Eliasz, P. (2005). Measuring the effects of
#' monetary policy: A factor-augmented vector autoregressive (FAVAR) approach.
#' \emph{The Quarterly Journal of Economics, 120}(1), 387--422.
#'
#' @export
create_favarmodel <- function(x, y, p = 2, n = 1, normalize_x = TRUE,
                              iterations = 20000, burnin = 2000) {

  # Input checks ----
  if (!"ts" %in% class(x)) {
    stop("Argument 'x' must be an object of class 'ts'.")
  }
  if (!"ts" %in% class(y)) {
    stop("Argument 'y' must be an object of class 'ts'.")
  }
  if (any(p < 0)) {
    stop("Argument 'p' must be at least 0.")
  }
  if (any(n < 1)) {
    stop("Argument 'n' must be at least 1.")
  }
  if (nrow(as.matrix(x)) != nrow(as.matrix(y))) {
    stop("Arguments 'x' and 'y' must cover the same periods.")
  }
  if (any(n >= NCOL(x))) {
    stop("Argument 'n' must be smaller than the number of columns of 'x': the ",
         "leading n series identify the factors and carry no free loading, so ",
         "there would be nothing left to estimate.")
  }

  # Data preparation ----
  if (normalize_x) {
    x <- scale(x)
  }

  m <- NCOL(x)
  n_obs <- NCOL(y)

  model <- NULL
  model$type <- "FAVAR"
  model$m <- m
  model$n <- 0
  model$n_obs <- n_obs
  model$p <- 0
  # The sampler add_posterior_coefficients dispatches on. Named rather than
  # derived at the point of use, so that the model object says which sampler
  # produced it.
  model$algorithm <- "FavarNormalWishart"
  model$iterations <- iterations
  model$burnin <- burnin

  result <- NULL
  for (j in n) {
    for (i in p) {
      model_i <- model
      model_i$n <- j
      model_i$p <- i

      result_i <- list("data" = list("x" = x, "y" = y),
                       "model" = model_i)

      class(result_i) <- append("favarmodel", class(result_i))

      result <- c(result, list(result_i))
    }
  }

  if (length(result) == 1) {
    result <- result[[1]]
  } else {
    class(result) <- append("modellist", class(result))
  }

  return(result)
}
