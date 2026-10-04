# expected_model_size() for the factor models.
#
# The size is predicted from the specification, so the invariant is that the
# prediction is what the samplers then store: for all four dynamic factor model
# samplers and the FAVAR, at one and at two factors, every element the
# prediction lists is in the posterior with exactly the predicted dimensions --
# the coefficients, the log likelihood and the forecasts -- and nothing in the
# posterior is missing from it. The warning built on the prediction belongs to
# bvartools and is tested there.

# The dimensions of every matrix in the posterior, as "element draws columns".
posterior_dims <- function(object) {
  post <- object[["posterior"]]
  result <- character(0)
  for (block in names(post)) {
    if (is.list(post[[block]])) {
      for (element in names(post[[block]])) {
        x <- post[[block]][[element]]
        result <- c(result, paste(paste0("posterior$", block, "$", element), NROW(x), NCOL(x)))
      }
    } else {
      result <- c(result, paste(paste0("posterior$", block), NROW(post[[block]]), NCOL(post[[block]])))
    }
  }
  sort(result)
}

predicted_dims <- function(size) {
  size <- size[size[["step"]] != "", ]
  sort(paste(size[["element"]], size[["draws"]], size[["columns"]]))
}

# Estimates the model, then compares every block with the prediction for the
# model as it stood before the step that added the block.
expect_size_matches <- function(object, label) {
  before <- expected_model_size(object)
  set.seed(1)
  fitted <- add_posterior_loglik(add_posterior_coefficients(object))
  fitted <- add_posterior_forecasts(add_forecast_input(fitted, n_ahead = 3))
  after <- expected_model_size(fitted)

  expect_identical(predicted_dims(after), posterior_dims(fitted), info = label)
  # Only the forecasts, whose horizon the forecast step sets, are unknown before.
  expect_identical(setdiff(predicted_dims(after), predicted_dims(before)),
                   grep("forecast", predicted_dims(after), value = TRUE), info = label)
}

test_that("the predicted blocks are what every factor model sampler stores", {
  for (error in c("gamma", "sv")) {
    for (tvp in c(FALSE, TRUE)) {
      for (n in 1:2) {
        sim <- sim_dfm(n = n, p = 2)
        object <- create_dfmodel(x = sim$x, p = 2, n = n, error = error, tvp = tvp,
                                 iterations = 7, burnin = 2)
        priors <- list(object)
        if (tvp) {
          priors[["lambda"]] <- tvp_prior()
          priors[["a"]] <- tvp_prior()
        }
        if (error == "sv") {
          priors[["u"]] <- sv_prior()
          priors[["v"]] <- sv_prior()
        }
        object <- add_initial_values(do.call(add_priors, priors))
        expect_size_matches(object, paste("DFM", error, tvp, n))
      }
    }
  }

  for (n in 1:2) {
    sim <- make_favar_sample()
    object <- create_favarmodel(x = sim$x, y = sim$yts, p = 2, n = n,
                                iterations = 7, burnin = 2)
    expect_size_matches(add_initial_values(add_priors(object)), paste("FAVAR", n))
  }
})

test_that("the size is eight bytes per number and a 'modelsize' to print", {
  object <- prepared_dfm_tvp(iterations = 7, burnin = 2)$object
  size <- expected_model_size(object)

  expect_s3_class(size, "modelsize")
  expect_equal(size[["bytes"]][-1], 8 * size[["draws"]][-1] * size[["columns"]][-1])
  expect_equal(size[["bytes"]][1],
               as.numeric(utils::object.size(object[names(object) != "posterior"])))
  expect_output(print(size), "posterior\\$lambda\\$coeffs")
})
