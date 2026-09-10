# The recursion the method uses is cheap because it propagates the impact
# vector rather than the matrices Psi_h. That is only correct if it agrees with
# the companion form, so the companion form is written out here and the two are
# compared draw by draw.
companion_response <- function(a_draw, q_draw, ns, p, impulse, h, order = NULL) {

  if (is.null(order)) {
    order <- seq_len(ns)
  }

  root <- matrix(0, ns, ns)
  root[order, order] <- t(chol(q_draw[order, order, drop = FALSE]))

  if (p == 0) {
    out <- matrix(0, ns, h + 1)
    out[, 1] <- root[, impulse]
    return(out)
  }

  a <- matrix(a_draw, ns, ns * p)
  companion <- matrix(0, ns * p, ns * p)
  companion[seq_len(ns), ] <- a
  if (p > 1) {
    companion[(ns + 1):(ns * p), seq_len(ns * (p - 1))] <- diag(ns * (p - 1))
  }

  out <- matrix(0, ns, h + 1)
  power <- diag(ns * p)
  for (i in 0:h) {
    out[, i + 1] <- (power %*% rbind(root, matrix(0, ns * (p - 1), ns)))[seq_len(ns), impulse]
    power <- power %*% companion
  }

  return(out)
}

q_of_draw <- function(object, i) {
  ns <- object$model$n + object$model$n_obs
  solve(matrix(object$posterior$v_sigma_inv$coeffs[i, ], ns, ns))
}

test_that("the response has the shape and the class a bvartools irf has", {

  est <- make_favar(iterations = 100, burnin = 50)

  ir <- irf(est$model, impulse = 1, response = 3, n_ahead = 8)

  expect_s3_class(ir, "bvarirf")
  expect_s3_class(ir, "ts")
  expect_equal(nrow(ir), 9)
  expect_equal(ncol(ir), 3)
  expect_equal(stats::start(ir)[1], 0)
  # The bands bracket the median, which is what makes the object readable.
  expect_true(all(ir[, 1] <= ir[, 2]))
  expect_true(all(ir[, 2] <= ir[, 3]))
})

test_that("the propagation agrees with the companion form", {

  est <- make_favar(iterations = 60, burnin = 40, p = 2)
  object <- est$model
  ns <- object$model$n + object$model$n_obs

  ir <- irf(object, impulse = 1, response = 1, n_ahead = 6, keep_draws = TRUE)
  ir <- as.matrix(ir)

  for (i in c(1, 7, nrow(ir))) {
    expected <- companion_response(object$posterior$a$coeffs[i, ],
                                   q_of_draw(object, i), ns, 2, 1, 6)
    expect_equal(as.numeric(ir[i, ]), expected[1, ], tolerance = 1e-10)
  }
})

test_that("the impact is the Cholesky factor and respects the ordering", {

  est <- make_favar(iterations = 60, burnin = 40)
  object <- est$model
  ns <- object$model$n + object$model$n_obs

  # Shock to the first element, ordered first: nothing else moves on impact,
  # and it moves by its own Cholesky diagonal.
  own <- as.matrix(irf(object, impulse = 1, response = 1, n_ahead = 0,
                       keep_draws = TRUE))
  other <- as.matrix(irf(object, impulse = 1, response = "obs", n_ahead = 0,
                         keep_draws = TRUE))

  root <- t(chol(q_of_draw(object, 1)))
  expect_equal(as.numeric(own[1, 1]), root[1, 1], tolerance = 1e-10)
  expect_true(all(abs(other) > 0))

  # Reversed, the observed block is ordered first and a shock to it is the one
  # that leaves the other unmoved on impact.
  rev_other <- as.matrix(irf(object, impulse = 1, response = "obs", n_ahead = 0,
                             order = c(2, 1), keep_draws = TRUE))
  expect_true(all(abs(rev_other) < 1e-12))
})

test_that("a shock scales the response linearly and cumulates on request", {

  est <- make_favar(iterations = 60, burnin = 40)
  object <- est$model

  one <- as.matrix(irf(object, impulse = 1, response = 2, n_ahead = 5,
                       keep_draws = TRUE))
  two <- as.matrix(irf(object, impulse = 1, response = 2, n_ahead = 5,
                       shock = -2, keep_draws = TRUE))
  expect_equal(two, -2 * one, tolerance = 1e-10)

  cumulated <- as.matrix(irf(object, impulse = 1, response = 2, n_ahead = 5,
                             cumulative = TRUE, keep_draws = TRUE))
  expect_equal(cumulated[3, ], cumsum(one[3, ]), tolerance = 1e-10)
})

test_that("a leading panel series responds exactly as its factor does", {

  # The identifying row of the loadings is a unit vector, so the panel series
  # that defines a factor and the factor itself cannot respond differently.
  est <- make_favar(iterations = 60, burnin = 40, n = 1)
  object <- est$model

  state <- as.matrix(irf(object, impulse = 1, response = "Series 1", n_ahead = 5,
                         keep_draws = TRUE))
  panel <- as.matrix(irf(object, impulse = 1, response = "Series 4", n_ahead = 5,
                         keep_draws = TRUE))

  lambda <- matrix(object$posterior$lambda$coeffs[1, ], object$model$m,
                   object$model$n + object$model$n_obs)
  expect_equal(lambda[1, ], c(1, 0), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(state, panel)))

  # And the panel response is the state's, carried through that series' row.
  s <- irf(object, impulse = 1, response = 1, n_ahead = 5, keep_draws = TRUE)
  y <- irf(object, impulse = 1, response = "obs", n_ahead = 5, keep_draws = TRUE)
  expect_equal(panel[1, ],
               lambda[4, 1] * as.matrix(s)[1, ] + lambda[4, 2] * as.matrix(y)[1, ],
               tolerance = 1e-10)
})

test_that("a transition of order zero responds on impact alone", {

  est <- make_favar(iterations = 40, burnin = 20, p = 0)
  object <- est$model

  expect_null(object$posterior$a)

  ir <- as.matrix(irf(object, impulse = 1, response = 1, n_ahead = 4,
                      keep_draws = TRUE))
  expect_true(all(abs(ir[, 1]) > 0))
  expect_true(all(abs(ir[, -1]) < 1e-12))
})

test_that("keep_draws returns one row per draw", {

  est <- make_favar(iterations = 60, burnin = 40)

  ir <- irf(est$model, impulse = 1, response = 2, n_ahead = 3,
            keep_draws = TRUE)

  expect_s3_class(ir, "mcmc")
  expect_equal(nrow(ir), 60)
  expect_equal(ncol(ir), 4)
})

test_that("the arguments are checked", {

  est <- make_favar(iterations = 40, burnin = 20)
  object <- est$model

  expect_error(irf(object, impulse = 1, response = 1),
               NA)
  expect_error(irf(object, impulse = "nowhere", response = 1),
               "not available as 'impulse'")
  expect_error(irf(object, impulse = 1, response = "nowhere"),
               "not available as 'response'")
  expect_error(irf(object, response = 1),
               "must be provided")
  # The impulse is a state element, so the panel's later series are out of range.
  expect_error(irf(object, impulse = 5, response = 1),
               "between 1 and 2")
  expect_error(irf(object, impulse = 1, response = 1, order = c(1, 1)),
               "permutation")
  expect_error(irf(object, impulse = 1, response = 1, ci = 1.5),
               "between 0 and 1")
  expect_error(irf(object, impulse = 1, response = 1, shock = c(1, 2)),
               "single number")

  no_posterior <- object
  no_posterior$posterior <- NULL
  expect_error(irf(no_posterior, impulse = 1, response = 1),
               "does not contain posterior draws")
})
