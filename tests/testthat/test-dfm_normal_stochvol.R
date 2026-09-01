# The stochastic volatility model, end to end through the R interface.
#
# What the numerics do is covered upstream, in BayesTS's
# test/unit_dfm_normal_stochvol.cpp. What is left for here is the translation:
# whether the right prior reaches the right error term, whether the two widths
# stay apart, and whether the wider posterior comes back with the shapes the
# documentation claims. The pair of widths is the thing worth testing on this
# side -- `u` is the number of observed series and `v` the number of factors, and
# a swap gives two well-formed lists.

test_that("create_dfmodel selects the stochastic volatility sampler", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1, p = 1)

  object <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv",
                           iterations = 20, burnin = 10)

  expect_s3_class(object, "dfmodel")
  expect_equal(object$model$error, "sv")
  expect_equal(object$model$algorithm, "DfmNormalStochvol")

  # And the gamma default is untouched.
  expect_equal(create_dfmodel(x = sim$x, p = 1, n = 1,
                              iterations = 20, burnin = 10)$model$algorithm,
               "DfmNormalGamma")

  expect_error(create_dfmodel(x = sim$x, p = 1, n = 1, error = "garch",
                              iterations = 20, burnin = 10),
               "Invalid specification")
})

test_that("add_priors builds both volatility blocks at their own width", {

  sim <- sim_dfm(tt = 40, m = 5, n = 2, p = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 2, error = "sv",
                           iterations = 20, burnin = 10)

  object <- add_priors(object, u = sv_prior(v_i = 0.5), v = sv_prior(v_i = 0.25))

  # u is M wide, v is N wide. Everything else about the two is the same.
  for (field in c("mu", "shape", "rate", "sigma", "offset")) {
    expect_equal(nrow(object$priors$u[[field]]), 5, info = field)
    expect_equal(nrow(object$priors$v[[field]]), 2, info = field)
  }
  expect_equal(dim(object$priors$u$v_inv), c(5, 5))
  expect_equal(dim(object$priors$v$v_inv), c(2, 2))

  # v_i arrives as a scalar precision and is stored as the matrix the core reads.
  expect_equal(diag(object$priors$u$v_inv), rep(0.5, 5))
  expect_equal(diag(object$priors$v$v_inv), rep(0.25, 2))

  # The coefficient blocks are the same as for the gamma model.
  expect_equal(dim(object$priors$lambda$vinv), rep(n_free_lambda(5, 2), 2))
  expect_equal(nrow(object$priors$a$mu), 2 * 2 * 1)
})

test_that("add_priors reports a gamma specification given to an sv model", {

  sim <- sim_dfm(tt = 40, m = 4, n = 1, p = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 1, error = "sv",
                           iterations = 20, burnin = 10)

  # The default arguments are the gamma ones, so this is what a caller who
  # forgets gets, and the message has to name what is missing.
  expect_error(add_priors(object), "'mu'")
  expect_error(add_priors(object, u = sv_prior(), v = list(shape = 5, rate = 4)),
               "'v' is missing")

  for (bad in list(list(v_i = 0), list(rate = 0), list(state_variance = 0),
                   list(offset = 0))) {
    spec <- utils::modifyList(sv_prior(), bad)
    expect_error(add_priors(object, u = spec, v = sv_prior()),
                 "must be larger than 0")
  }
})

test_that("add_initial_values adds a flat log-volatility path per error term", {

  sim <- sim_dfm(tt = 40, m = 5, n = 2, p = 1)
  object <- create_dfmodel(x = sim$x, p = 1, n = 2, error = "sv",
                           iterations = 20, burnin = 10)
  object <- add_priors(object, u = sv_prior(), v = sv_prior())

  set.seed(7)
  object <- add_initial_values(object)

  expect_equal(dim(object$initial$u_h), c(40, 5))
  expect_equal(dim(object$initial$v_h), c(40, 2))
  expect_equal(dim(object$initial$u_h_init), c(5, 1))
  expect_equal(dim(object$initial$v_h_init), c(2, 1))

  # Flat, and starting exactly at the state before the sample.
  expect_equal(object$initial$u_h[1, ], object$initial$u_h[40, ])
  expect_equal(object$initial$u_h[1, ], as.numeric(object$initial$u_h_init))
  expect_equal(object$initial$v_h[1, ], as.numeric(object$initial$v_h_init))

  # The gamma model's precisions are not there, and the sv model's paths are not
  # there for the gamma model.
  expect_null(object$initial$uinv)
  expect_null(object$initial$vinv)

  gamma_object <- add_initial_values(add_priors(
    create_dfmodel(x = sim$x, p = 1, n = 2, iterations = 20, burnin = 10)))
  expect_null(gamma_object$initial$u_h)
  expect_equal(dim(gamma_object$initial$uinv), c(5, 5))
})

test_that("add_posterior_coefficients returns a volatility path per draw", {

  m <- 4
  n <- 2
  p <- 2
  tt <- 40
  iterations <- 20

  prep <- prepared_dfm_sv(iterations = iterations, burnin = 10,
                          tt = tt, m = m, n = n, p = p)

  set.seed(1)
  object <- add_posterior_coefficients(prep$object)

  expect_s3_class(object, "dfmodel")
  expect_false(isTRUE(object$error))
  expect_named(object$posterior,
               c("lambda", "factors", "a", "u_sigma_inv", "v_sigma_inv"))

  for (i in names(object$posterior)) {
    expect_s3_class(object$posterior[[i]]$coeffs, "mcmc")
    expect_equal(nrow(object$posterior[[i]]$coeffs), iterations)
  }

  # The coefficient blocks are the same shapes DfmNormalGamma returns.
  expect_equal(ncol(object$posterior$lambda$coeffs), m * n)
  expect_equal(ncol(object$posterior$factors$coeffs), tt * n)
  expect_equal(ncol(object$posterior$a$coeffs), n * n * p)

  # The error blocks are not: a whole path per draw rather than one number per
  # series, so each widens by a factor of tt.
  expect_equal(ncol(object$posterior$u_sigma_inv$coeffs), m * tt)
  expect_equal(ncol(object$posterior$v_sigma_inv$coeffs), n * tt)

  # Precisions, so positive, and finite throughout.
  expect_true(all(object$posterior$u_sigma_inv$coeffs > 0))
  expect_true(all(object$posterior$v_sigma_inv$coeffs > 0))
  expect_true(all(is.finite(object$posterior$factors$coeffs)))
  expect_true(all(is.finite(object$posterior$lambda$coeffs)))

  # And the volatility actually moves. A path that came back constant would mean
  # the per-period precision never reached the output, which every shape check
  # above would still pass.
  uvar <- idiosyncratic_variance(object, m, tt)
  expect_true(stats::sd(colMeans(uvar)) > 1e-8)

  # The identifying block is still fixed in every draw.
  for (i in seq_len(nrow(object$posterior$lambda$coeffs))) {
    draw <- matrix(object$posterior$lambda$coeffs[i, ], m, n)
    block <- draw[seq_len(n), , drop = FALSE]
    expect_equal(diag(block), rep(1, n))
    expect_true(all(block[upper.tri(block)] == 0))
  }
})

test_that("the sampler draws from R's RNG, so set.seed reproduces the draws", {

  # The same guard the gamma model's tests carry, for the second binding: losing
  # the RcppArmadillo wiring compiles, links, runs and silently stops honouring
  # set.seed().
  prep <- prepared_dfm_sv(iterations = 20, burnin = 10, tt = 40, m = 4, n = 2, p = 1)

  set.seed(11)
  first <- add_posterior_coefficients(prep$object)$posterior

  set.seed(11)
  second <- add_posterior_coefficients(prep$object)$posterior

  expect_equal(first, second)

  set.seed(12)
  third <- add_posterior_coefficients(prep$object)$posterior

  expect_false(isTRUE(all.equal(as.numeric(first$u_sigma_inv$coeffs),
                                as.numeric(third$u_sigma_inv$coeffs))))
})

test_that("the estimated volatility follows a break in the data", {

  # The point of the model, at the coarsest resolution that is still an assertion
  # about the numbers: a sample whose idiosyncratic noise quadruples half way
  # through has to come back with a higher variance in the second half.
  # DfmNormalGamma cannot express this at all -- one variance per series for the
  # whole sample -- so it is also what distinguishes the two.
  m <- 4
  n <- 1
  tt <- 200

  sim <- sim_dfm(tt = tt, m = m, n = n, p = 1)
  set.seed(99)
  half <- tt %/% 2
  noise <- matrix(stats::rnorm(tt * m, sd = 0.3), tt, m)
  noise[(half + 1):tt, ] <- noise[(half + 1):tt, ] * 4
  x <- stats::ts(sim$f %*% t(sim$lambda) + noise,
                 start = stats::start(sim$x), frequency = 4)
  colnames(x) <- colnames(sim$x)

  object <- create_dfmodel(x = x, p = 1, n = n, normalize_x = FALSE,
                           error = "sv", iterations = 150, burnin = 150)
  object <- add_priors(object, u = sv_prior(), v = sv_prior())
  set.seed(5)
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  expect_false(isTRUE(object$error))

  uvar <- idiosyncratic_variance(object, m, tt)
  early <- mean(uvar[, seq_len(half)])
  late <- mean(uvar[, (half + 1):tt])

  # A quadrupling of the standard deviation is a sixteenfold variance; asserting
  # a factor of three leaves room for the smoothing the random walk imposes and
  # for a short chain, while a flat path would fail it outright.
  expect_gt(late / early, 3)
})

test_that("forecasts and the log-likelihood dispatch to the sv sampler", {

  m <- 4
  n <- 2
  tt <- 40
  iterations <- 20
  n_ahead <- 4

  prep <- prepared_dfm_sv(iterations = iterations, burnin = 10,
                          tt = tt, m = m, n = n, p = 1)

  set.seed(3)
  object <- add_posterior_coefficients(prep$object)
  object <- add_posterior_forecasts(object, n_ahead = n_ahead)

  expect_s3_class(object$posterior$forecast, "mcmc")
  expect_equal(dim(object$posterior$forecast), c(iterations, n_ahead * m))
  expect_true(all(is.finite(object$posterior$forecast)))

  object <- add_posterior_loglik(object)
  expect_equal(dim(object$posterior$loglik), c(iterations, tt))
  expect_true(all(is.finite(object$posterior$loglik)))

  # Both refuse an object without the factor path, which neither can recompute.
  no_factors <- object
  no_factors$posterior$factors <- NULL
  expect_error(add_posterior_forecasts(no_factors, n_ahead = 2), "factors")
  expect_error(add_posterior_loglik(no_factors), "factors")

  # And an unknown algorithm is reported rather than silently taken for the
  # gamma one, which is what the dispatch added to these two is for.
  wrong <- object
  wrong$model$algorithm <- "DfmSomethingElse"
  expect_error(add_posterior_forecasts(wrong, n_ahead = 2), "not supported")
  expect_error(add_posterior_loglik(wrong), "not supported")
})

test_that("a mismatched volatility width is refused by the sampler", {

  # The failure this is here for: the two blocks have different widths, so a
  # caller who fills one from the other gets a well-formed list of the wrong
  # length. add_posterior_coefficients traps the error and returns error = TRUE,
  # so the assertion is on that rather than on a condition.
  prep <- prepared_dfm_sv(iterations = 20, burnin = 10, tt = 40, m = 4, n = 2, p = 1)

  object <- prep$object
  object$initial$v_h <- object$initial$u_h # 40 x 4, where 40 x 2 is wanted

  expect_true(suppressWarnings(add_posterior_coefficients(object))$error)
})
