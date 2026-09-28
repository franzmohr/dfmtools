# Sign and zero restrictions on the state equation of a factor augmented VAR.
#
# The algorithm is bvartools' and is tested there. What is checked here is the
# translation: that the draws handed to it are the FAVAR's transition and
# state covariance, that every identified draw satisfies what was asked of it
# once irf() reads it back, that the posterior is resampled as a whole, and
# that what cannot be restricted is refused.

sz_fit <- function() {
  cached <- getOption("dfmtools.test.sz_fit")
  if (!is.null(cached)) return(cached)
  fit <- make_favar(iterations = 300, burnin = 100, p = 1)
  options(dfmtools.test.sz_fit = fit)
  fit
}

sz_state <- function(model) {
  c(colnames(model[["data"]][["x"]])[1], colnames(model[["data"]][["y"]]))
}

test_that("the identified responses satisfy the restrictions", {
  model <- sz_fit()$model
  state <- sz_state(model)
  # The shock named after the factor raises it and leaves the observed
  # variable unmoved on impact; the observed variable's own shock raises it.
  # The zero goes on the shock drawn first: the last column of the rotation
  # has no dimension left to spare for one.
  restrictions <- data.frame(impulse = c(state[1], state[1], state[2]),
                             response = c(state[1], state[2], state[2]),
                             sign = c(1, 0, 1))
  set.seed(1)
  identified <- add_sign_zero_restrictions(model, restrictions, max_tries = 50)

  record <- identified[["model"]][["sign_zero_restrictions"]]
  expect_gt(record[["accepted"]], 0)
  expect_identical(record[["candidates"]], nrow(model[["posterior"]][["a"]][["coeffs"]]))

  k <- ncol(model[["data"]][["x"]])
  obs_on_factor <- irf(identified, impulse = 1, response = k + 1, type = "sign",
                       n_ahead = 3, keep_draws = TRUE)
  obs_on_obs <- irf(identified, impulse = 2, response = k + 1, type = "sign",
                    n_ahead = 3, keep_draws = TRUE)
  factor_on_factor <- irf(identified, impulse = 1, response = 1, type = "sign",
                          n_ahead = 3, keep_draws = TRUE)
  expect_equal(max(abs(obs_on_factor[, 1])), 0, tolerance = 1e-8)
  expect_true(all(obs_on_obs[, 1] > 0))
  expect_true(all(factor_on_factor[, 1] > 0))
})

test_that("the posterior is resampled as a whole", {
  model <- sz_fit()$model
  state <- sz_state(model)
  set.seed(2)
  identified <- add_sign_zero_restrictions(
    model, data.frame(impulse = state[1], response = state, sign = c(1, 0)), draws = 40)

  for (block in c("lambda", "factors", "a", "u_sigma_inv", "v_sigma_inv", "q")) {
    draws <- identified[["posterior"]][[block]][["coeffs"]]
    expect_identical(nrow(draws), 40L, info = block)
    expect_s3_class(draws, "mcmc")
  }
  # Each row of every block is one original draw: the transition kept with a
  # rotation is the one it was found for.
  original <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])
  kept <- as.matrix(identified[["posterior"]][["a"]][["coeffs"]])
  v_original <- as.matrix(model[["posterior"]][["v_sigma_inv"]][["coeffs"]])
  v_kept <- as.matrix(identified[["posterior"]][["v_sigma_inv"]][["coeffs"]])
  row <- match(kept[1, 1], original[, 1])
  expect_equal(v_kept[1, ], v_original[row, ])
})

test_that("the variance decomposition reads the rotations and adds up", {
  model <- sz_fit()$model
  state <- sz_state(model)
  set.seed(3)
  identified <- add_sign_zero_restrictions(
    model, data.frame(impulse = state[1], response = state, sign = c(1, 0)))
  k <- ncol(model[["data"]][["x"]])

  vd <- fevd(identified, response = k + 1, n_ahead = 4, type = "sign")
  expect_equal(unname(rowSums(vd)), rep(1, 5), tolerance = 1e-8)
  # A panel series has an idiosyncratic share as well, and the three add up.
  panel <- fevd(identified, response = 3, n_ahead = 2, type = "sign")
  expect_equal(unname(rowSums(panel)), rep(1, 3), tolerance = 1e-8)
})

test_that("signs alone are accepted", {
  model <- sz_fit()$model
  state <- sz_state(model)
  set.seed(4)
  identified <- add_sign_zero_restrictions(
    model, data.frame(impulse = state, response = state, sign = c(1, 1)), max_tries = 20)
  expect_gt(identified[["model"]][["sign_zero_restrictions"]][["accepted"]], 0)
})

test_that("what cannot be restricted is refused with a reason", {
  model <- sz_fit()$model
  state <- sz_state(model)
  panel <- colnames(model[["data"]][["x"]])
  expect_error(add_sign_zero_restrictions(
    model, data.frame(impulse = state[2], response = panel[3], sign = 1)),
    "responds through its loadings")
  expect_error(irf(model, impulse = 1, response = 1, type = "sign"),
               "add_sign_zero_restrictions")
  expect_error(fevd(model, response = 1, type = "sign"), "add_sign_zero_restrictions")
  expect_error(irf(model, impulse = 1, response = 1, type = "sign", order = 2:1),
               "only meaningful")
})

test_that("an exactly identified rotation is the Cholesky factor it must be", {
  # With two state elements, a zero from the factor's shock to the observed
  # variable makes the impact matrix upper triangular, and the two signs fix
  # its columns: the identification is the Cholesky one with the observed
  # variable ordered first, draw by draw. A wrong convention for turning the
  # rotation into an impact matrix -- a transpose, the other triangle -- fails
  # this where the sign checks above would still pass.
  model <- sz_fit()$model
  state <- sz_state(model)
  restrictions <- data.frame(impulse = c(state[1], state[1], state[2]),
                             response = c(state[1], state[2], state[2]),
                             sign = c(1, 0, 1))
  set.seed(5)
  identified <- add_sign_zero_restrictions(model, restrictions, max_tries = 50)
  k <- ncol(model[["data"]][["x"]])

  for (impulse in 1:2) {
    for (response in c(1, 3, k + 1)) {
      sign <- irf(identified, impulse = impulse, response = response, type = "sign",
                  n_ahead = 4, keep_draws = TRUE)
      cholesky <- irf(identified, impulse = impulse, response = response, type = "oir",
                      order = 2:1, n_ahead = 4, keep_draws = TRUE)
      expect_equal(unclass(sign), unclass(cholesky), tolerance = 1e-8, ignore_attr = TRUE)
    }
  }
})
