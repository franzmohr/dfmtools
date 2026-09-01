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
