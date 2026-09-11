test_that("add_initial_values produces one starting value per drawn block", {

  sim <- sim_dfm(tt = 40, m = 4, n = 2)
  object <- add_initial_values(add_priors(create_dfmodel(x = sim$x, p = 2, n = 2)))

  m <- 4
  n <- 2
  p <- 2

  expect_named(object$initial, c("lambda", "uinv", "vinv", "a"))

  # The free loadings only, in the column-by-column order this side stores them
  # in -- the binding permutes them to the row-major order the core reads.
  expect_equal(dim(object$initial$lambda), c(n_free_lambda(m, n), 1L))
  expect_equal(dim(object$initial$a), c(n * n * p, 1L))

  # Precisions, and diagonal: only the diagonal of either covariance matrix is
  # estimated under the gamma specification.
  expect_equal(dim(object$initial$uinv), c(m, m))
  expect_equal(dim(object$initial$vinv), c(n, n))
  expect_true(all(object$initial$uinv[upper.tri(object$initial$uinv)] == 0))
  expect_true(all(object$initial$vinv[upper.tri(object$initial$vinv)] == 0))
  expect_true(all(diag(object$initial$uinv) > 0))
  expect_true(all(diag(object$initial$vinv) > 0))
})

test_that("a transition of order zero gets no starting value for a", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- add_initial_values(add_priors(create_dfmodel(x = sim$x, p = 0, n = 1)))

  expect_null(object$initial$a)
  expect_false(is.null(object$initial$lambda))
})

test_that("starting values are drawn at the scale of the prior, not its inverse", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  # Every prior here is a precision, so a draw from it has standard deviation
  # 1/sqrt(vinv) -- which is what backsolve(chol(vinv), z) gives. The other
  # arrangement of the same two objects, chol(vinv) %*% z, gives sqrt(vinv), and
  # this file held both at once until they were reconciled.
  #
  # The two agree at vinv = 1 and nowhere else, and they move in opposite
  # directions as the prior tightens, so the check runs at a diffuse precision
  # and at a tight one. A tight prior is the case that matters: there the wrong
  # form throws the starting values far outside the prior rather than merely
  # bunching them inside it.
  scale_of <- function(vinv, reps = 1000) {
    object <- add_priors(create_dfmodel(x = sim$x, p = 1, n = 1),
                         lambda = list(vinv = vinv), a = list(vinv = vinv))
    draws <- replicate(reps, {
      # method = "prior", because the loadings are the one block the default
      # does not draw -- see the mirror mode tests below.
      initial <- add_initial_values(object, method = "prior")$initial
      c(initial$lambda[1], initial$a[1])
    })
    apply(draws, 1, stats::sd)
  }

  set.seed(20)
  expect_equal(scale_of(0.01), c(10, 10), tolerance = 0.1)
  expect_equal(scale_of(100), c(0.1, 0.1), tolerance = 0.1)
})

test_that("starting values are drawn from R's RNG, so set.seed fixes them", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- add_priors(create_dfmodel(x = sim$x, p = 1, n = 1))

  set.seed(7)
  first <- add_initial_values(object)$initial

  set.seed(7)
  second <- add_initial_values(object)$initial

  expect_equal(first, second)
})

# The loadings are the one block whose starting value decides which mode the
# sampler ends up in rather than how long it takes to get there. Only the
# product lambda %*% f_t is identified and the restriction fixes the leading
# N x N block of lambda, so a start whose free loadings have the wrong sign is a
# coherent model -- the mirror, in which the factor is the negative of the common
# component and the identifying series is treated as noise -- and the chain stays
# in it. Hence method = "pca": the principal components estimate is what the data
# say, and the rotation onto the restriction makes it independent of the
# arbitrary sign svd() returns each column with.
#
# The packed vector holds the free elements column by column, each column from
# the diagonal down, which is what this rebuilds to compare against.
rebuild_lambda <- function(packed, m, n) {
  lambda <- diag(1, m, n)
  at <- 1
  for (j in seq_len(n)) {
    for (i in seq_len(m - j) + j) {
      lambda[i, j] <- packed[at]
      at <- at + 1
    }
  }
  lambda
}

pca_lambda <- function(x, n) {
  v <- svd(unclass(x), nu = 0, nv = n)$v
  v %*% solve(v[seq_len(n), , drop = FALSE])
}

test_that("the loadings start where the principal components put them", {

  for (spec in list(c(m = 4, n = 1), c(m = 5, n = 2), c(m = 6, n = 3), c(m = 3, n = 3))) {

    m <- spec[["m"]]
    n <- spec[["n"]]
    sim <- sim_dfm(tt = 80, m = m, n = n, p = 1)
    object <- add_initial_values(add_priors(create_dfmodel(x = sim$x, p = 1, n = n)))

    expect_equal(dim(object$initial$lambda), c(n_free_lambda(m, n), 1L), info = paste(m, n))
    expect_equal(rebuild_lambda(object$initial$lambda, m, n),
                 pca_lambda(object$data$x, n), info = paste(m, n))
  }
})

test_that("the loading start does not depend on the seed, and the rest still does", {

  sim <- sim_dfm(tt = 80, m = 5, n = 2)
  object <- add_priors(create_dfmodel(x = sim$x, p = 1, n = 2))

  initial_at <- function(seed, ...) {
    set.seed(seed)
    add_initial_values(object, ...)$initial
  }

  expect_equal(initial_at(1)$lambda, initial_at(2)$lambda)
  expect_false(isTRUE(all.equal(initial_at(1)$a, initial_at(2)$a)))

  # Drawing them is still available, and still a draw.
  expect_false(isTRUE(all.equal(initial_at(1, method = "prior")$lambda,
                                initial_at(2, method = "prior")$lambda)))
})

test_that("a time varying loading path starts flat at the same estimate", {

  sim <- sim_dfm(tt = 80, m = 5, n = 2, p = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 2, tvp = TRUE)
  object <- add_priors(object, lambda = tvp_prior(), a = tvp_prior())
  object <- add_initial_values(object)

  k <- n_free_lambda(5, 2)
  expect_equal(dim(object$initial$lambda), c(k, sim$tt))
  expect_equal(rebuild_lambda(object$initial$lambda[, 1], 5, 2),
               pca_lambda(object$data$x, 2))

  # Flat, and at the pre-sample state, as every drifting block here starts.
  expect_equal(object$initial$lambda[, sim$tt], object$initial$lambda[, 1])
  expect_equal(object$initial$lambda_init[, 1], object$initial$lambda[, 1])
})

test_that("a panel the estimate cannot be rotated for falls back to the prior", {

  # Two identical series in the identifying positions leave the leading block of
  # the component loadings singular, so there is no rotation onto the
  # restriction and nothing data-based to start from.
  sim <- sim_dfm(tt = 80, m = 5, n = 2)
  x <- sim$x
  x[, 2] <- x[, 1]

  object <- add_priors(create_dfmodel(x = x, p = 1, n = 2))

  set.seed(3)
  expect_warning(initial <- add_initial_values(object)$initial,
                 "too little of the common variation")
  expect_true(all(is.finite(initial$lambda)))

  # It is the prior draw, which is what the same seed gives when asked for it.
  set.seed(3)
  expect_equal(initial$lambda, add_initial_values(object, method = "prior")$initial$lambda)
})

test_that("only 'pca' and 'prior' are methods", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)
  object <- add_priors(create_dfmodel(x = sim$x, p = 1, n = 1))

  expect_error(add_initial_values(object, method = "ols"), "'pca' or 'prior'")
})

test_that("the sampler does not mirror the factor from the default start", {

  # Six series whose loadings are known, including two negative ones, over a
  # sample long enough for the mode to be decided by where the chain starts.
  # Every one of these seeds lands on a factor that moves with the panel; the
  # same sweep under method = "prior" mirrors on some of them, which is the
  # behaviour this default exists to remove.
  sim <- sim_dfm_known(tt = 200)

  correlation <- vapply(1:6, function(seed) {
    set.seed(seed)
    object <- create_dfmodel(x = sim$x, p = 1, n = 1, iterations = 300, burnin = 200)
    object <- add_priors(object)
    object <- add_initial_values(object)
    object <- add_posterior_coefficients(object)
    stats::cor(colMeans(object$posterior$factors$coeffs), sim$x[, 1])
  }, numeric(1))

  expect_true(all(correlation > 0.9))
})

# One observed series and one factor is the only specification whose loading
# matrix is entirely fixed: lambda is 1 x 1, the identifying restriction puts a
# one in it, and (2M - N - 1)N/2 is zero. What is left is x_t = f_t + u_t with an
# AR(p) factor, which is an unobserved-components model and which the sampler
# estimates like any other -- so the empty block is carried through rather than
# rejected. Before this was handled, chol() met the 0 x 0 prior precision and
# reported a zero-dimensional matrix instead of anything about the model.
single_series <- function(tt = 60) {
  set.seed(11)
  f <- numeric(tt + 20)
  for (i in 2:(tt + 20)) f[i] <- 0.7 * f[i - 1] + stats::rnorm(1)
  f <- f[-(1:20)]
  stats::ts(matrix(f + stats::rnorm(tt, sd = 0.5), ncol = 1,
                   dimnames = list(NULL, "a")),
            start = c(1980, 1), frequency = 4)
}

test_that("one series and one factor leaves an empty loading block, not an error", {

  x <- single_series()

  object <- add_priors(create_dfmodel(x = x, p = 1, n = 1,
                                      iterations = 20, burnin = 10))
  expect_equal(n_free_lambda(1, 1), 0)
  expect_equal(dim(object$priors$lambda$vinv), c(0L, 0L))

  object <- add_initial_values(object)

  expect_equal(dim(object$initial$lambda), c(0L, 1L))
  expect_equal(dim(object$initial$a), c(1L, 1L))
  expect_equal(dim(object$initial$uinv), c(1L, 1L))
  expect_equal(dim(object$initial$vinv), c(1L, 1L))

  # And the sampler runs it, which is why the specification is supported rather
  # than turned away. The whole of lambda comes back, as it does for every other
  # model here, and every draw of it is the fixed one.
  object <- add_posterior_coefficients(object)
  expect_equal(dim(object$posterior$lambda$coeffs), c(20L, 1L))
  expect_true(all(object$posterior$lambda$coeffs == 1))
  expect_equal(dim(object$posterior$a$coeffs), c(20L, 1L))
})

test_that("the empty loading block is empty for time varying coefficients too", {

  x <- single_series()
  tt <- nrow(x)

  object <- create_dfmodel(x = x, p = 1, n = 1, tvp = TRUE,
                           iterations = 20, burnin = 10)
  object <- add_priors(object, lambda = tvp_prior(), a = tvp_prior())

  # The state equation of the loadings has nothing to describe either.
  expect_length(object$priors$lambda$shape, 0)
  expect_length(object$priors$lambda$rate, 0)

  object <- add_initial_values(object)

  expect_equal(dim(object$initial$lambda), c(0L, tt))
  expect_equal(dim(object$initial$lambda_init), c(0L, 1L))
  expect_equal(dim(object$initial$lambda_sigma_inv), c(0L, 0L))

  # The transition still drifts, and is the only block that does.
  expect_equal(dim(object$initial$a), c(1L, tt))

  object <- add_posterior_coefficients(object)
  expect_equal(dim(object$posterior$lambda$coeffs), c(20L, tt))
  expect_true(all(object$posterior$lambda$coeffs == 1))
})
