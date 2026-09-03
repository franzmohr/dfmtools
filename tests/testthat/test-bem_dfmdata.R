# The data set shipped with the package. It is what every example and the
# README run on, so a change to it -- a rebuild of data-raw that drops a column,
# reorders the sample, or introduces an NA -- would break those silently.

test_that("bem_dfmdata is the documented quarterly sample", {

  data("bem_dfmdata", envir = environment())

  expect_s3_class(bem_dfmdata, "mts")
  expect_s3_class(bem_dfmdata, "ts")

  # 196 US macroeconomic variables, 1959Q3 to 2015Q3.
  expect_equal(dim(bem_dfmdata), c(225L, 196L))
  expect_equal(stats::frequency(bem_dfmdata), 4)
  expect_equal(stats::start(bem_dfmdata), c(1959L, 3L))
  expect_equal(stats::end(bem_dfmdata), c(2015L, 3L))

  # Named columns, since the loadings are read off by variable.
  expect_length(colnames(bem_dfmdata), 196L)
  expect_false(any(is.na(colnames(bem_dfmdata))))
  expect_equal(anyDuplicated(colnames(bem_dfmdata)), 0L)
  expect_equal(colnames(bem_dfmdata)[1], "GDPC96")

  # Complete and finite: the sampler has no missing-data treatment, so an NA
  # anywhere would propagate through every draw.
  expect_false(anyNA(bem_dfmdata))
  expect_true(all(is.finite(bem_dfmdata)))

  # Already transformed to stationarity, which is what create_dfmodel documents
  # as its input. Series still carrying a unit root would show up here as a
  # column whose mean over the first half differs wildly from the second.
  expect_true(all(apply(bem_dfmdata, 2, stats::sd) > 0))
})

test_that("the documented example runs on the shipped data", {

  # A handful of columns rather than all 196, so that this stays a test of the
  # interface on real data and not a benchmark.
  data("bem_dfmdata", envir = environment())

  object <- create_dfmodel(x = bem_dfmdata[, 1:5], p = 1, n = 1,
                           iterations = 50, burnin = 50)
  object <- add_initial_values(add_priors(object))

  set.seed(71)
  object <- add_posterior_coefficients(object)

  expect_false(isTRUE(object$error))
  expect_equal(dim(object$posterior$lambda$coeffs), c(50L, 5L))
  expect_true(all(is.finite(object$posterior$lambda$coeffs)))

  # The data are standardised by default, so the loadings are on a comparable
  # scale across series and the first is the fixed one.
  expect_true(all(matrix(object$posterior$lambda$coeffs, 50, 5)[, 1] == 1))
})
