test_that("create_dfmodel returns a dfmodel carrying the specification", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  object <- create_dfmodel(x = sim$x, p = 2, n = 1,
                           iterations = 100, burnin = 50)

  expect_s3_class(object, "dfmodel")
  expect_named(object, c("data", "model"))

  expect_equal(object$model$type, "DFM")
  expect_equal(object$model$m, 4L)
  expect_equal(object$model$n, 1)
  expect_equal(object$model$p, 2)
  expect_equal(object$model$error, "gamma")
  expect_equal(object$model$iterations, 100)
  expect_equal(object$model$burnin, 50)

  # The sampler add_posterior_coefficients dispatches on is chosen here, not
  # guessed there.
  expect_equal(object$model$algorithm, "DfmNormalGamma")
})

test_that("normalize_x controls whether the columns are standardised", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  normalised <- create_dfmodel(x = sim$x, p = 1, n = 1)$data$x
  expect_equal(dim(normalised), c(40L, 4L))
  expect_equal(colMeans(normalised), setNames(rep(0, 4), colnames(sim$x)))
  expect_equal(apply(normalised, 2, stats::sd), setNames(rep(1, 4), colnames(sim$x)))

  raw <- create_dfmodel(x = sim$x, p = 1, n = 1, normalize_x = FALSE)$data$x
  expect_equal(as.numeric(raw), as.numeric(sim$x))
})

test_that("integer vectors for p and n produce every combination", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  object <- create_dfmodel(x = sim$x, p = 1:2, n = 1:3)

  expect_s3_class(object, "modellist")
  expect_length(object, 6)
  expect_true(all(vapply(object, inherits, logical(1), "dfmodel")))

  specs <- t(vapply(object, function(z) c(z$model$p, z$model$n), numeric(2)))
  expect_equal(specs[order(specs[, 2], specs[, 1]), ],
               as.matrix(expand.grid(p = 1:2, n = 1:3))[, c("p", "n")],
               ignore_attr = TRUE)
})

test_that("a single specification is returned unwrapped", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  object <- create_dfmodel(x = sim$x, p = 1, n = 1)

  expect_s3_class(object, "dfmodel")
  expect_false(inherits(object, "modellist"))
})

test_that("as many factors as series is allowed, more than that is not", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  # N = M is a real specification: the identifying block is then the whole
  # leading square of lambda, and the N(N - 1)/2 elements below its diagonal are
  # still freely estimated. The sampler runs it.
  expect_s3_class(create_dfmodel(x = sim$x, p = 1, n = 4), "dfmodel")

  # N above M is not, because the unit lower triangle that pins the rotation and
  # scale of the factors needs one row of lambda per factor to be pinned with.
  # Caught here rather than left to add_priors(), where the count of free
  # loadings it implies is merely a wrong number -- and, far enough above M, a
  # negative one that surfaces as diag() refusing a negative dimension.
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 5),
               "must not exceed the number of observed series, which is 4")
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 20),
               "cannot have more factors than series")

  # And the same for a vector of n, where one entry reaching past M is enough.
  expect_s3_class(create_dfmodel(x = sim$x, p = 1, n = 1:4), "modellist")
  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1:5),
               "must not exceed the number of observed series")
})

test_that("a univariate ts is turned into a one-column matrix and named", {

  set.seed(11)
  uni <- stats::ts(stats::rnorm(40), start = c(1980, 1), frequency = 4)

  object <- create_dfmodel(x = uni, p = 1, n = 1, iterations = 20, burnin = 10)

  x <- object$data$x

  # A univariate ts arrives as a vector, and everything downstream counts on a
  # matrix with a column to name.
  expect_equal(dim(x), c(40L, 1L))
  expect_equal(colnames(x), "y")
  expect_equal(object$model$m, 1L)

  # The conversion must not cost the series its time index, which is what
  # as.matrix() would drop if it were not carried over.
  expect_equal(stats::tsp(x), stats::tsp(uni))
})

test_that("a univariate ts keeps its values and index unnormalised", {

  set.seed(12)
  uni <- stats::ts(stats::rnorm(40), start = c(1980, 1), frequency = 4)

  x <- create_dfmodel(x = uni, p = 1, n = 1, normalize_x = FALSE)$data$x

  expect_equal(as.numeric(x), as.numeric(uni))
  expect_equal(stats::tsp(x), stats::tsp(uni))
  expect_equal(colnames(x), "y")
})

test_that("a ts matrix without column names is given default ones", {

  set.seed(13)
  mv <- stats::ts(matrix(stats::rnorm(120), 40, 3), start = c(1980, 1), frequency = 4)
  # ts() names the columns "Series 1" and so on when the matrix it is given has
  # none, so they have to be taken off again to get one that really carries no
  # names.
  dimnames(mv) <- NULL

  object <- create_dfmodel(x = mv, p = 1, n = 1)

  x <- object$data$x

  # Named rather than left empty, so that nothing downstream has to work with
  # columns it cannot refer to.
  expect_equal(colnames(x), c("y1", "y2", "y3"))
  expect_equal(dim(x), c(40L, 3L))
  expect_equal(object$model$m, 3L)
  expect_equal(stats::tsp(x), stats::tsp(mv))
})

test_that("a ts matrix with row names but no column names is named too", {

  set.seed(14)
  mv <- stats::ts(matrix(stats::rnorm(80), 40, 2), start = c(1980, 1), frequency = 4)
  dimnames(mv) <- list(paste0("obs", seq_len(40)), NULL)

  x <- create_dfmodel(x = mv, p = 1, n = 1)$data$x

  expect_equal(colnames(x), c("y1", "y2"))
})

test_that("column names that are there are left alone", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  x <- create_dfmodel(x = sim$x, p = 1, n = 1)$data$x

  expect_equal(colnames(x), colnames(sim$x))
})

test_that("an unnamed ts matrix carries through the whole preparation", {

  set.seed(15)
  mv <- stats::ts(matrix(stats::rnorm(80), 40, 2), start = c(1980, 1), frequency = 4)
  dimnames(mv) <- NULL

  object <- create_dfmodel(x = mv, p = 1, n = 1, iterations = 20, burnin = 10)
  object <- add_priors(object)
  object <- add_initial_values(object)

  set.seed(16)
  object <- add_posterior_coefficients(object)

  expect_false(isTRUE(object$error))
  expect_equal(ncol(object$posterior$lambda$coeffs), 2)
  expect_equal(ncol(object$posterior$u_sigma_inv$coeffs), 2)

})

test_that("a univariate ts carries through to a posterior", {

  set.seed(17)
  uni <- stats::ts(stats::rnorm(60), start = c(1980, 1), frequency = 4)

  object <- create_dfmodel(x = uni, p = 1, n = 1, iterations = 20, burnin = 10)
  object <- add_priors(object)
  object <- add_initial_values(object)

  set.seed(18)
  object <- add_posterior_coefficients(object)

  expect_false(isTRUE(object$error))

  # One series with one factor leaves no freely estimated loading, so this
  # reaches the sampler at all only because the empty block is carried through --
  # single_series() in test-add_initial_values.R covers that specification given
  # as a one-column matrix. What is tested here is the vector, which is the form
  # a univariate ts actually arrives in.
  expect_equal(dim(object$posterior$lambda$coeffs), c(20L, 1L))
  expect_true(all(object$posterior$lambda$coeffs == 1))
  expect_equal(ncol(object$posterior$factors$coeffs), 60)
})

test_that("create_dfmodel rejects invalid input", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1)

  expect_error(create_dfmodel(x = as.data.frame(sim$x)), "must be an object of class 'ts'")
  expect_error(create_dfmodel(x = matrix(as.numeric(sim$x), nrow(sim$x))),
               "must be an object of class 'ts'")
  expect_error(create_dfmodel(x = sim$x, p = -1), "'p' must be at least 0")
  expect_error(create_dfmodel(x = sim$x, n = 0), "'n' must be at least 1")
  expect_error(create_dfmodel(x = sim$x, error = "wishart"), "Invalid specification")
  expect_error(create_dfmodel(x = sim$x, error = 1), "must be of class 'character'")
})
