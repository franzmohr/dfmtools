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
