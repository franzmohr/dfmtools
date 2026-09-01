# A small dynamic factor model, simulated rather than taken from bem_dfmdata:
# 40 quarters of 4 variables keeps the sampler in these tests to a fraction of a
# second, and the true loadings and transition are known if a test ever needs
# them.
sim_dfm <- function(tt = 40, m = 4, n = 1, p = 1, seed = 42) {

  set.seed(seed)

  # Transition: a stationary diagonal A_1, further lags at zero.
  a <- diag(0.5, n)

  # Burn 20 periods in so the path does not start at the origin.
  f <- matrix(0, tt + 20, n)
  for (i in 2:(tt + 20)) {
    f[i, ] <- a %*% f[i - 1, ] + stats::rnorm(n, sd = 0.5)
  }
  f <- f[-(1:20), , drop = FALSE]

  # Loadings with the identifying block the sampler fixes: unit lower triangular
  # in the leading n rows.
  lambda <- matrix(stats::runif(m * n, 0.5, 1.5), m, n)
  lead <- diag(n)
  lead[lower.tri(lead)] <- stats::runif(n * (n - 1) / 2)
  lambda[seq_len(n), ] <- lead

  x <- f %*% t(lambda) + matrix(stats::rnorm(tt * m, sd = 0.5), tt, m)
  colnames(x) <- paste0("x", seq_len(m))

  list(x = stats::ts(x, start = c(1980, 1), frequency = 4),
       f = f, lambda = lambda, a = a, p = p, n = n, m = m, tt = tt)
}

# create_dfmodel -> add_priors -> add_initial_values, the state the sampler needs.
prepared_dfm <- function(iterations = 20, burnin = 10, ...) {

  sim <- sim_dfm(...)

  object <- create_dfmodel(x = sim$x, p = sim$p, n = sim$n,
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object)
  object <- add_initial_values(object)

  list(object = object, sim = sim)
}

# Number of freely estimated loadings, i.e. the M x N matrix less the fixed ones
# and zeros of the identifying block.
n_free_lambda <- function(m, n) (2 * m - n - 1) * n / 2

# A stochastic volatility prior for one of the two error terms. add_priors() has
# the gamma specification as its default, so an sv model has to be given this;
# the widths are filled in from the model, so one list serves both u and v.
sv_prior <- function(mu = 0, v_i = 1, shape = 3, rate = 0.2,
                     state_variance = 0.05, offset = 1e-4) {
  list(mu = mu, v_i = v_i, shape = shape, rate = rate,
       state_variance = state_variance, offset = offset)
}

# prepared_dfm() for error = "sv".
prepared_dfm_sv <- function(iterations = 20, burnin = 10, ...) {

  sim <- sim_dfm(...)

  object <- create_dfmodel(x = sim$x, p = sim$p, n = sim$n, error = "sv",
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object, u = sv_prior(), v = sv_prior())
  object <- add_initial_values(object)

  list(object = object, sim = sim)
}

# The idiosyncratic variances a posterior implies, M x T: the stored object is a
# precision path, m values per period, periods along a row.
idiosyncratic_variance <- function(object, m, tt) {
  matrix(1 / colMeans(object$posterior$u_sigma_inv$coeffs), m, tt)
}
