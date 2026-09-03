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

# The state equation one coefficient block of a time varying model needs, on top
# of the `vinv` add_priors() already takes for it. There is no default for the
# pair, so a tvp model has to be given this for both lambda and a.
tvp_prior <- function(vinv = 0.01, shape = 3, rate = 0.01) {
  list(vinv = vinv, shape = shape, rate = rate)
}

# prepared_dfm() for tvp = TRUE.
prepared_dfm_tvp <- function(iterations = 20, burnin = 10, ...) {

  sim <- sim_dfm(...)

  object <- create_dfmodel(x = sim$x, p = sim$p, n = sim$n, tvp = TRUE,
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object, lambda = tvp_prior(), a = tvp_prior())
  object <- add_initial_values(object)

  list(object = object, sim = sim)
}

# prepared_dfm() for tvp = TRUE and error = "sv", which is the two specifications
# at once: a state equation for each coefficient block and a volatility block for
# each error term.
prepared_dfm_tvp_sv <- function(iterations = 20, burnin = 10, ...) {

  sim <- sim_dfm(...)

  object <- create_dfmodel(x = sim$x, p = sim$p, n = sim$n, error = "sv", tvp = TRUE,
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object, lambda = tvp_prior(), a = tvp_prior(),
                       u = sv_prior(), v = sv_prior())
  object <- add_initial_values(object)

  list(object = object, sim = sim)
}

# The loading path a posterior implies, (M N) x T: the stored object is the whole
# M x N matrix per period, periods along a row. At one factor the rows are the M
# series in order.
loading_path <- function(object, m, n, tt) {
  matrix(colMeans(object$posterior$lambda$coeffs), m * n, tt)
}

# A sample whose loadings really move: one factor, and every series' loading on a
# ramp from `from` to `to` over the sample. What a constant-loading model cannot
# find, and what makes a recovery test of this model worth running.
#
# The first loading is 1 at both ends because the identifying restriction fixes
# it there -- see sim_dfm_known(), which says why at length. The others move by
# enough that a sampler that ignored the period index would fail on the level as
# well as on the direction.
sim_dfm_drifting <- function(tt = 400, u_sd = 0.4, u_sd_to = u_sd, a = 0.6, seed = 7,
                             from = c(1, 1.5, -0.8, 2.0, 0.5),
                             to   = c(1, 0.3, -0.8, 0.6, 1.5)) {

  set.seed(seed)

  m <- length(from)

  # Burn 50 periods in so the path does not start at the origin.
  f <- numeric(tt + 50)
  for (i in 2:(tt + 50)) {
    f[i] <- a * f[i - 1] + stats::rnorm(1)
  }
  f <- f[-(1:50)]

  share <- (seq_len(tt) - 1) / (tt - 1)
  lambda <- outer(1 - share, from) + outer(share, to) # tt x m

  # The idiosyncratic scale ramps too when `u_sd_to` differs from `u_sd`, which
  # is what a model carrying both drifts has to tell apart from the loadings
  # moving. With the default the two are equal and the scale is constant.
  u_sd_t <- (1 - share) * u_sd + share * u_sd_to

  x <- lambda * f + matrix(stats::rnorm(tt * m), tt, m) * u_sd_t
  colnames(x) <- paste0("x", seq_len(m))

  list(x = stats::ts(x, start = c(1950, 1), frequency = 4),
       f = f, lambda = lambda, from = from, to = to,
       u_sd = u_sd, u_sd_to = u_sd_to, u_sd_t = u_sd_t,
       n = 1, m = m, tt = tt)
}

# The idiosyncratic variances a posterior implies, M x T: the stored object is a
# precision path, m values per period, periods along a row.
idiosyncratic_variance <- function(object, m, tt) {
  matrix(1 / colMeans(object$posterior$u_sigma_inv$coeffs), m, tt)
}

# A sample with parameters that are known exactly, for the tests that assert on
# the numbers rather than on the shapes. sim_dfm() draws its loadings at random
# and gives every series the same noise, neither of which a recovery test can
# check itself against.
#
# The first loading is 1 because the identifying restriction fixes it there:
# only the product lambda %*% f_t is identified, so a factor scaled by c and
# loadings scaled by 1/c is the same model. Pinning the true value to the value
# the sampler holds fixed is what makes the remaining loadings comparable to
# truth directly, without rescaling the draw first.
sim_dfm_known <- function(tt = 300,
                          lambda = matrix(c(1, 1.5, -0.8, 2, 0.5, -1.2), 6, 1),
                          u_sd = c(0.3, 0.5, 0.4, 0.6, 0.3, 0.5),
                          a = 0.6, seed = 123) {

  set.seed(seed)

  m <- nrow(lambda)
  n <- ncol(lambda)
  a <- diag(a, n)

  # Burn 50 periods in so the path does not start at the origin.
  f <- matrix(0, tt + 50, n)
  for (i in 2:(tt + 50)) {
    f[i, ] <- a %*% f[i - 1, ] + stats::rnorm(n, sd = 1)
  }
  f <- f[-(1:50), , drop = FALSE]

  # A different noise scale per series, so that the idiosyncratic variances are
  # distinguishable from one another and a sampler that returned a common one
  # would fail.
  x <- f %*% t(lambda) + matrix(stats::rnorm(tt * m, sd = rep(u_sd, each = tt)), tt, m)
  colnames(x) <- paste0("x", seq_len(m))

  list(x = stats::ts(x, start = c(1950, 1), frequency = 4),
       f = f, lambda = lambda, a = a, u_sd = u_sd,
       n = n, m = m, tt = tt)
}
