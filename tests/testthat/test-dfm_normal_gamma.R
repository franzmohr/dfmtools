# The gamma model, end to end through the R interface.
#
# test-add_posterior_coefficients.R covers the shapes this sampler returns and
# the wiring that returns them. What is left for here is whether the numbers
# mean anything: a sampler that has lost a transpose, indexed a loading column
# by row, or dropped the prior scale still returns an mcmc object of exactly the
# documented width. The sample is small enough (300 quarters, 6 series, one
# factor) that a chain long enough to have converged runs in about a second.

test_that("the sampler recovers loadings, variances and the transition", {

  sim <- sim_dfm_known(tt = 300)

  # normalize_x = FALSE, because the recovery is asserted against the loadings
  # the data were generated with, and standardising the columns rescales them.
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, normalize_x = FALSE,
                           iterations = 1000, burnin = 1000)
  object <- add_initial_values(add_priors(object))

  set.seed(7)
  object <- add_posterior_coefficients(object)

  expect_false(isTRUE(object$error))

  m <- sim$m
  n <- sim$n
  tt <- sim$tt

  # The loadings. The first is the fixed one, so it is an assertion about the
  # identifying restriction rather than about the draw; the other five are
  # estimated, and they span both signs and a fourfold range of magnitude, so a
  # sampler that returned the prior mean or a permutation of the column would
  # miss them.
  lambda_hat <- matrix(colMeans(object$posterior$lambda$coeffs), m, n)
  expect_equal(lambda_hat[1, 1], 1)
  expect_equal(as.numeric(lambda_hat), as.numeric(sim$lambda), tolerance = 0.25)

  # The idiosyncratic variances, which differ by series in the sample. The
  # stored block is a precision, so this inverts it. The estimates sit above the
  # truth -- the factor path is itself estimated, so part of the common
  # variation lands in the idiosyncratic term -- by around 0.04 on variances of
  # 0.09 to 0.36. That is a large relative difference on the smallest of them,
  # so the bound is absolute rather than the relative tolerance expect_equal()
  # would apply.
  u_var <- 1 / colMeans(object$posterior$u_sigma_inv$coeffs)
  expect_true(all(abs(u_var - sim$u_sd^2) < 0.15))

  # A sampler returning one variance for every series would pass the check above
  # only by accident; the spread across series has to survive too.
  expect_gt(stats::cor(u_var, sim$u_sd^2), 0.9)

  # The transition coefficient of the factor.
  expect_equal(as.numeric(colMeans(object$posterior$a$coeffs)),
               as.numeric(diag(sim$a)), tolerance = 0.2)

  # And the factor path itself, which is drawn rather than computed. Its scale
  # is pinned by the fixed loading, so this compares to the simulated path
  # directly.
  f_hat <- t(matrix(colMeans(object$posterior$factors$coeffs), n, tt))
  expect_gt(stats::cor(f_hat[, 1], sim$f[, 1]), 0.95)
})

test_that("the fitted measurement equation reproduces the sample", {

  # The residuals of the posterior mean fit have to be small relative to the
  # data and centred at zero. This is what ties the two halves of a draw
  # together: lambda and the factor path can each be individually plausible and
  # still not multiply back out to x if the reshape of either is wrong.
  sim <- sim_dfm_known(tt = 300)

  object <- create_dfmodel(x = sim$x, p = 1, n = 1, normalize_x = FALSE,
                           iterations = 500, burnin = 500)
  object <- add_initial_values(add_priors(object))

  set.seed(8)
  object <- add_posterior_coefficients(object)

  lambda_hat <- matrix(colMeans(object$posterior$lambda$coeffs), sim$m, sim$n)
  f_hat <- t(matrix(colMeans(object$posterior$factors$coeffs), sim$n, sim$tt))

  fitted <- f_hat %*% t(lambda_hat)
  residuals <- sim$x - fitted

  expect_equal(colMeans(residuals), rep(0, sim$m), tolerance = 0.1,
               ignore_attr = TRUE)

  # The common component explains most of each series: the idiosyncratic noise
  # was built at a fraction of the loading, so a residual variance anywhere near
  # the variance of the series itself would mean the factor is not being fitted.
  expect_true(all(apply(residuals, 2, stats::var) < 0.5 * apply(sim$x, 2, stats::var)))
})

test_that("a transition of order zero runs through the whole pipeline", {

  # p = 0 is accepted by create_dfmodel and leaves the transition equation with
  # nothing to estimate. add_priors and add_initial_values are already known to
  # skip their `a` blocks; what is untested is that the sampler and the two
  # functions downstream of it do the same rather than reading an absent block.
  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  object <- create_dfmodel(x = sim$x, p = 0, n = 1, iterations = 20, burnin = 10)
  object <- add_initial_values(add_priors(object))

  set.seed(51)
  object <- add_posterior_coefficients(object)

  expect_false(isTRUE(object$error))
  expect_null(object$posterior$a$coeffs)

  # Everything else is drawn at its usual width.
  expect_equal(ncol(object$posterior$lambda$coeffs), 4L)
  expect_equal(ncol(object$posterior$factors$coeffs), 40L)
  expect_equal(ncol(object$posterior$u_sigma_inv$coeffs), 4L)
  expect_equal(ncol(object$posterior$v_sigma_inv$coeffs), 1L)

  # A factor with no dynamics still has a forecast: it is white noise around
  # zero, and the measurement equation still maps it onto M variables.
  object <- add_posterior_forecasts(object, n_ahead = 3)
  expect_equal(dim(object$posterior$forecast), c(20L, 3L * 4L))
  expect_true(all(is.finite(object$posterior$forecast)))

  object <- add_posterior_loglik(object)
  expect_equal(dim(object$posterior$loglik), c(20L, 40L))
  expect_true(all(is.finite(object$posterior$loglik)))
})
