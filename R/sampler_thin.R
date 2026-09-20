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


# An argument that has to be a vector of whole numbers at or above `minimum`.
# `p` and `n` of create_dfmodel() and create_favarmodel(), which index a model
# rather than measure anything.
#
# The `any(p < 0)` this replaces let a character through: "a" < 0 compares as
# strings and is FALSE, so a typed argument travelled on and failed several
# steps later, inside a matrix multiplication, in terms of a quantity the caller
# never named.
.check_model_integer <- function(value, name, minimum) {
  if (!is.numeric(value) || length(value) == 0 || anyNA(value) ||
      any(!is.finite(value)) || any(value != round(value))) {
    stop("Argument '", name, "' must be a vector of whole numbers.", call. = FALSE)
  }
  if (any(value < minimum)) {
    stop("Argument '", name, "' must be at least ", minimum, ".", call. = FALSE)
  }
  return(invisible(value))
}

# The length of a chain, validated where the caller names it rather than where
# the sampler trips over it. The core refuses a non-positive `iterations` and a
# negative `burnin` too, but only once add_posterior_coefficients() runs, which
# is a priors and starting values round trip after the mistake was made -- and
# there the message arrives through the `error = TRUE` path rather than as a
# stop. `thin` has been checked here all along; these two now are as well.
.check_sampler_length <- function(value, name, minimum) {
  if (!is.numeric(value) || length(value) != 1 || is.na(value) ||
      !is.finite(value) || value != round(value) || value < minimum) {
    stop("Argument '", name, "' must be a single whole number of at least ",
         minimum, ".", call. = FALSE)
  }
  return(invisible(value))
}
