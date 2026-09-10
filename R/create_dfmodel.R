#' Create a Dynamic Factor Model
#'
#' Produces the input for the estimation of a dynamic factor model (DFM).
#'
#' @param x a time-series object of stationary endogenous variables.
#' @param p an integer vector of the lag order of the measurement equation. See 'Details'.
#' @param n an integer vector of the number of factors, which must not exceed the number
#' of columns of \code{x}. See 'Details'.
#' @param normalize_x logical indicating whether each column of \code{x} should
#' be normalized using \code{scale}. Defaults to \code{TRUE}.
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"gamma"}. See 'Details'.
#' @param tvp logical indicating whether the coefficients of the model are time varying.
#' Defaults to \code{FALSE}. See 'Details'.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 20000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#'
#' @details The function produces the variable matrices of dynamic factor
#' models (DFM) with measurement equation
#' \deqn{x_t = \lambda f_t + u_t,}
#' where
#' \eqn{x_t} is an \eqn{M \times 1} vector of observed variables,
#' \eqn{f_t} is an \eqn{N \times 1} vector of unobserved factors and
#' \eqn{\lambda} is the corresponding \eqn{M \times N} matrix of factor loadings.
#' \eqn{u_t} is an \eqn{M \times 1} error term with \eqn{u_t \sim N(0, U)}.
#'
#' The transition equation is
#' \deqn{f_t = \sum_{i=1}^{p} A_i f_{t - i} + v_t,}
#' where
#' \eqn{A_i} is an \eqn{N \times N} coefficient matrix and
#' \eqn{v_t} is an \eqn{N \times 1} error term with \eqn{v_t \sim N(0, V)}.
#'
#' If integer vectors are provided as arguments \code{p} or \code{n}, the function will
#' produce a distinct model for all possible combinations of those specifications.
#'
#' Only the product \eqn{\lambda f_t} is identified, so the leading \eqn{N \times N} block of
#' \eqn{\lambda} is fixed unit lower triangular -- ones on the diagonal, zeros above -- which pins
#' both the rotation and the scale of the factors. That leaves \eqn{N(2M - N - 1)/2} freely
#' estimated loadings and requires \eqn{N \le M}, since the block needs one row of \eqn{\lambda}
#' per factor to be fixed with; a larger \code{n} is rejected here rather than left to the sampler.
#' \eqn{N = M} is allowed, and the \eqn{N(N - 1)/2} elements below the diagonal of the leading
#' block remain free. The one specification with no free loading at all is \eqn{M = N = 1}, where
#' \eqn{\lambda} is the single fixed one and the model reduces to \eqn{x_t = f_t + u_t} with an
#' AR(\eqn{p}) factor; that is estimated like any other, with an empty loading block throughout.
#'
#' Argument \code{error} specifies the structure of the covariance matrix of
#' the error term and how it is estimated. Possible specifications are:
#' \itemize{
#'  \item{\code{"gamma"}: Only the diagonal elements of the covariance matrix are estimated using a gamma prior.
#' Off-diagonal elements are not estimated and set to zero.}
#'  \item{\code{"sv"}: Only the diagonal elements of the covariance matrix are estimated, and they move
#' with time: the log-volatility of every error term follows a random walk of its own, estimated with the
#' mixture approximation of Omori et al. (2007). This applies to both error terms, \eqn{u_t} of the
#' measurement equation and \eqn{v_t} of the transition equation, and the two do different things.
#' Volatility in \eqn{u_t} reweights the series that identify the factors, so a series that was noisy
#' early and quiet later stops contributing on the same terms throughout. Volatility in \eqn{v_t} is the
#' common component's own, and it is what keeps the \eqn{M} idiosyncratic variances from jointly absorbing
#' a shock that every series felt at once. Off-diagonal elements are not estimated and set to zero.}
#' }
#'
#' If \code{tvp = TRUE}, the coefficients of the model follow random walks and are estimated as
#' whole paths with the simulation smoother of Durbin and Koopman (2002). \emph{Both} coefficient
#' blocks move: every freely estimated element of the loading matrix \eqn{\lambda_t} of the
#' measurement equation, and every element of the coefficient matrices \eqn{A_{i, t}} of the
#' transition equation. A loading is a series' exposure to the common component, and that it held
#' over the whole sample is the assumption a factor model makes most often and defends least: a
#' constant-loading model has nowhere to put a change in exposure except the idiosyncratic
#' variance, which then carries it as noise the series is credited with throughout. The leading
#' \eqn{N \times N} block of \eqn{\lambda_t} does not move, because it is not estimated -- it is
#' the normalisation that identifies the factors, and letting it drift would let their rotation and
#' scale wander over the sample.
#'
#' The two arguments combine, and the four algorithms they select are
#' \code{DfmNormalGamma}, \code{DfmNormalStochvol}, \code{DfmTvpGamma} and
#' \code{DfmTvpStochvol}. Carrying both at once is not redundant: a series whose loading fell
#' looks like a series whose idiosyncratic variance rose, and a period of common turbulence looks
#' like a transition that changed, so a model with only one of the two has to explain the other
#' with what it has.
#'
#' @return An object of class \code{'dfmodel'}, which contains the following elements:
#' \item{data}{A list of data objects, which can be used for posterior simulation. Element
#' \code{X} is a time-series object of normalised observable variables, i.e. each column has
#' zero mean and unity variance.}
#' \item{model}{A list of model specifications.}
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model data
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1,
#'                         iterations = 5000, burnin = 1000)
#'
#' # The same with stochastic volatility in both error terms
#' model_sv <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1, error = "sv",
#'                            iterations = 5000, burnin = 1000)
#'
#' # And with time varying loadings and transition coefficients
#' model_tvp <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1, tvp = TRUE,
#'                             iterations = 5000, burnin = 1000)
#'
#' # Both at once
#' model_tvp_sv <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1, error = "sv", tvp = TRUE,
#'                                iterations = 5000, burnin = 1000)
#'
#' @references
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). \emph{Bayesian Econometric Methods}
#' (2nd ed.). Cambridge: University Press.
#'
#' Del Negro, M., & Otrok, C. (2008). Dynamic factor models with time-varying parameters: measuring
#' changes in international business cycles. \emph{Federal Reserve Bank of New York Staff Report}
#' No. 326.
#'
#' Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation smoother for state space
#' time series analysis. \emph{Biometrika 89}(3), 603--615.
#'
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#'
#' Omori, Y., Chib, S., Shephard, N., & Nakajima, J. (2007). Stochastic volatility with leverage.
#' Fast and efficient likelihood inference. \emph{Journal of Econometrics 140}(2), 425--449.
#'
#' @export
create_dfmodel <- function(x, p = 2, n = 1, normalize_x = TRUE, error = "gamma", tvp = FALSE, iterations = 20000, burnin = 2000) {
  
  
  # Input checks ----
  if (!"ts" %in% class(x)) {
    stop("Argument 'x' must be an object of class 'ts'.")
  }
  
  if (any(p < 0)) {
    stop("Argument 'p' must be at least 0.")
  }
  
  if (any(n < 1)) {
    stop("Argument 'n' must be at least 1.")
  }
  
  if ("character" %in% class(error)) {
    if (!error %in% c("gamma", "sv")) {
      stop("Invalid specification of argument 'error'.")
    }
  } else {
    stop("Argument 'error' must be of class 'character'.")
  }
  
  if (!"logical" %in% class(tvp)) {
    stop("Argument 'tvp' must be of class 'logical'.")
  }
  
  
  # Data preparation ----
  
  if (is.null(dimnames(x))) {
    tsp_temp <- stats::tsp(x)
    data <- stats::ts(as.matrix(x), class = c("mts", "ts", "matrix"))
    stats::tsp(x) <- tsp_temp
    dimnames(x)[[2]] <- "y"
  }
  
  # Normalise every column of x
  if (normalize_x) {
    x <- scale(x)
  }
  
  data_name <- dimnames(x)[[2]]
  m <- NCOL(x)
  tt <- nrow(x)
  p_max <- max(p)

  # More factors than series is not a model. Only the product lambda f_t is
  # identified, so the leading N x N block of lambda is fixed unit lower
  # triangular, and there has to be one row per factor to fix it with. Rejected
  # here rather than left to the sampler, because the count of free loadings
  # that add_priors() derives from an N above M is a wrong number rather than an
  # impossible one -- at N far enough above M it turns negative, and what a
  # caller then sees is diag() complaining about a negative dimension.
  #
  # N equal to M is allowed and does run: the identifying block is then the
  # whole of the leading square, and N(N - 1)/2 loadings are still free below it.
  if (any(n > m)) {
    stop("Argument 'n' must not exceed the number of observed series, which is ",
         m, ". A dynamic factor model cannot have more factors than series.")
  }
  
  model <- NULL
  model$type <- "DFM"
  model$m <- m
  model$n <- 0
  model$p <- 0
  model$error <- error
  model$tvp <- tvp
  # The sampler add_posterior_coefficients dispatches on. Named rather than
  # derived from `error` and `tvp` at the point of use, so that the model object
  # says which sampler produced it.
  model$algorithm <- if (tvp) {
    switch(error,
           "gamma" = "DfmTvpGamma",
           "sv" = "DfmTvpStochvol")
  } else {
    switch(error,
           "gamma" = "DfmNormalGamma",
           "sv" = "DfmNormalStochvol")
  }
  model$iterations <- iterations
  model$burnin <- burnin
  
  
  
  result <- NULL
  for (j in n) {
    for (i in p) {
      model_i <- model
      model_i$n <- j
      model_i$p <- i
      
      result_i <- list("data" = list("x" = x),
                       "model" = model_i)
      
      class(result_i) <- append("dfmodel", class(result_i))
      
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
