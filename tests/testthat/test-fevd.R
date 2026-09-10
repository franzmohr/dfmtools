# The generalised responses and the variance decomposition. Both are checked
# against identities that hold exactly rather than against stored numbers: a
# generalised response coincides with an orthogonalised one for the element
# ordered first, and an orthogonalised decomposition sums to one.

test_that("a generalised response equals the orthogonalised one for the first element", {

  # Q e_1 / sqrt(q_11) is the first column of the Cholesky factor, so the two
  # identifications cannot disagree about a shock to the element ordered first.
  est <- make_favar(iterations = 80, burnin = 40)
  object <- est$model

  oir <- as.matrix(irf(object, impulse = 1, response = 2, n_ahead = 6,
                       keep_draws = TRUE))
  gir <- as.matrix(irf(object, impulse = 1, response = 2, n_ahead = 6,
                       type = "gir", keep_draws = TRUE))
  expect_equal(gir, oir, tolerance = 1e-10)

  # For the element ordered last they differ, the shock no longer being the
  # first Cholesky column.
  oir2 <- as.matrix(irf(object, impulse = 2, response = 2, n_ahead = 6,
                        keep_draws = TRUE))
  gir2 <- as.matrix(irf(object, impulse = 2, response = 2, n_ahead = 6,
                        type = "gir", keep_draws = TRUE))
  expect_false(isTRUE(all.equal(gir2, oir2)))

  # And a generalised response does not depend on the ordering at all, which is
  # the whole point of it: reversing the state leaves it where it was.
  rev_oir <- as.matrix(irf(object, impulse = 2, response = 2, n_ahead = 6,
                           order = c(2, 1), keep_draws = TRUE))
  expect_equal(rev_oir, gir2, tolerance = 1e-10)
})

test_that("a forecast error response ignores the covariance entirely", {

  est <- make_favar(iterations = 60, burnin = 40)
  object <- est$model

  # A unit shock to element 1 leaves element 2 unmoved on impact, whatever Q
  # says, because no orthogonalisation is applied.
  feir <- as.matrix(irf(object, impulse = 1, response = "obs", n_ahead = 0,
                        type = "feir", keep_draws = TRUE))
  expect_true(all(abs(feir) < 1e-12))

  own <- as.matrix(irf(object, impulse = 1, response = 1, n_ahead = 0,
                       type = "feir", keep_draws = TRUE))
  expect_true(all(abs(own - 1) < 1e-12))
})

test_that("an orthogonalised decomposition sums to one at every horizon", {

  est <- make_favar(iterations = 60, burnin = 40)
  object <- est$model

  # An observed variable has no idiosyncratic term, so the state shocks alone
  # account for the whole forecast error.
  vd <- fevd(object, response = "obs", n_ahead = 8)
  expect_s3_class(vd, "bvarfevd")
  expect_equal(ncol(vd), 2)
  expect_equal(colnames(vd), c("Series 1", "obs"))
  expect_equal(as.numeric(rowSums(vd)), rep(1, 9), tolerance = 1e-10)

  # A panel series carries one more column, and the total still sums to one.
  vd_panel <- fevd(object, response = 4, n_ahead = 8)
  expect_equal(ncol(vd_panel), 3)
  expect_equal(colnames(vd_panel)[3], "idiosyncratic")
  expect_equal(as.numeric(rowSums(vd_panel)), rep(1, 9), tolerance = 1e-10)
  expect_true(all(vd_panel >= 0))
})

test_that("the idiosyncratic share is what the common component leaves over", {

  est <- make_favar(iterations = 60, burnin = 40)
  object <- est$model

  # Series 1 identifies the factor, so its loading is a unit vector and its
  # forecast error is the factor's plus its own noise.
  vd <- fevd(object, response = 1, n_ahead = 6)
  expect_true(all(vd[, "idiosyncratic"] > 0))
  # The idiosyncratic share falls with the horizon: the common component
  # accumulates, the white noise does not.
  expect_true(all(diff(vd[, "idiosyncratic"]) < 0))
})

test_that("a generalised decomposition overlaps unless normalised", {

  est <- make_favar(iterations = 60, burnin = 40)
  object <- est$model

  plain <- fevd(object, response = 4, n_ahead = 6, type = "gir")
  # The generalised shocks are not orthogonal, so the shares do not sum to one.
  expect_false(isTRUE(all.equal(as.numeric(rowSums(plain)), rep(1, 7))))

  normalised <- fevd(object, response = 4, n_ahead = 6, type = "gir",
                     normalise_gir = TRUE)
  expect_equal(as.numeric(rowSums(normalised)), rep(1, 7), tolerance = 1e-10)
  # Normalising leaves the idiosyncratic share alone, that one being exact.
  expect_equal(normalised[, "idiosyncratic"], plain[, "idiosyncratic"],
               tolerance = 1e-12)
})

test_that("the decomposition of a p = 0 model is its impact alone", {

  est <- make_favar(iterations = 40, burnin = 20, p = 0)

  vd <- fevd(est$model, response = "obs", n_ahead = 3)
  # With no transition, every horizon looks like the first.
  expect_equal(as.numeric(vd[1, ]), as.numeric(vd[4, ]), tolerance = 1e-10)
  expect_equal(as.numeric(rowSums(vd)), rep(1, 4), tolerance = 1e-10)
})

test_that("the arguments are checked", {

  est <- make_favar(iterations = 40, burnin = 20)
  object <- est$model

  expect_error(irf(object, impulse = 1, response = 1, type = "nonsense"),
               "must be one of")
  expect_error(irf(object, impulse = 1, response = 1, type = "gir",
                   order = c(2, 1)),
               "only meaningful for type")
  expect_error(fevd(object, response = 1, type = "feir"),
               "not defined")
  expect_error(fevd(object, response = 1, type = "gir", order = c(2, 1)),
               "only meaningful for type")
  expect_error(fevd(object, response = "nowhere"),
               "not available as 'response'")
  expect_error(fevd(object), "must be provided")

  no_posterior <- object
  no_posterior$posterior <- NULL
  expect_error(fevd(no_posterior, response = 1),
               "does not contain posterior draws")
})
