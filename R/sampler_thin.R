# The thinning interval a model's sampler keeps one draw in, validated. NULL for
# one, which is how a model that keeps every draw says so.
.check_sampler_thin <- function(thin) {
  if (!is.numeric(thin) || length(thin) != 1 || is.na(thin) ||
      thin < 1 || thin != round(thin)) {
    stop("Argument 'thin' must be a single positive integer.")
  }
  if (thin == 1) {
    return(NULL)
  }
  return(as.integer(thin))
}


# The sampler's draws of one block as an mcmc object whose labels say which draws
# were kept. BayesTS keeps the last of every `thin` draws after the burn-in, so the
# kept ones are iterations thin, 2 thin, ... after the burn-in. Without thinning
# this is exactly what coda::as.mcmc() returns, so an unthinned model is unchanged
# down to the type of its labels. The same helper as in bvartools, which this
# package cannot reach, since it is not exported there.
.mcmc_draws <- function(model, draws) {
  thin <- model[["thin"]]
  if (is.null(thin) || thin == 1) {
    return(coda::as.mcmc(draws))
  }
  thin <- as.numeric(thin)
  return(coda::mcmc(draws, start = thin, end = NROW(draws) * thin, thin = thin))
}
