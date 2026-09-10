# A factor augmented VAR whose factor is known, so recovery can be checked
# rather than only shape.
#
# The panel is built so that the first series *is* the factor plus noise, which
# is what the identification asserts: the leading n x n block of the factor
# loadings is the identity and the observed columns of those rows are zero. That
# is the restriction to test, because it is the one that differs from a dynamic
# factor model's and the one a wrong choice would leave running and plausible.
make_favar_sample <- function(tt = 160, n_x = 8, noise = 0.3, seed = 42) {
  set.seed(seed)
  f <- as.numeric(stats::filter(stats::rnorm(tt), 0.7, method = "recursive"))
  y <- as.numeric(stats::filter(stats::rnorm(tt), 0.5, method = "recursive"))

  # Series 1 loads one on the factor and zero on the observed block: the
  # identifying row.
  load_f <- c(1, seq(0.8, 1.3, length.out = n_x - 1))
  load_y <- c(0, seq(0.2, 0.6, length.out = n_x - 1))

  x <- outer(f, load_f) + outer(y, load_y) +
    matrix(stats::rnorm(tt * n_x, sd = noise), tt, n_x)

  list(f = f,
       y = y,
       load_f = load_f,
       load_y = load_y,
       x = stats::ts(x, start = c(1990, 1), frequency = 4),
       yts = stats::ts(matrix(y, dimnames = list(NULL, "obs")),
                       start = c(1990, 1), frequency = 4))
}

# normalize_x = FALSE so that the estimated loadings are on the same scale as
# the ones the sample was built with. With the default normalisation each panel
# column is divided by its own standard deviation, so a loading comes back scaled
# by that -- correct, and not comparable to the truth without undoing it.
make_favar <- function(iterations = 600, burnin = 200, p = 2, n = 1) {
  sim <- make_favar_sample()
  model <- create_favarmodel(x = sim$x, y = sim$yts, p = p, n = n,
                             normalize_x = FALSE,
                             iterations = iterations, burnin = burnin)
  model <- add_priors(model)
  model <- add_initial_values(model)
  list(sim = sim, model = add_posterior_coefficients(model))
}

