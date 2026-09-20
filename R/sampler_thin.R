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

# A single whole number at or above `minimum`: the length of a chain, or the
# horizon of a forecast or a response. Validated where the caller names it
# rather than where whatever consumes it trips over it. The core refuses a non-positive `iterations` and a
# negative `burnin` too, but only once add_posterior_coefficients() runs, which
# is a priors and starting values round trip after the mistake was made -- and
# there the message arrives through the `error = TRUE` path rather than as a
# stop. `thin` has been checked here all along; these two now are as well.
.check_whole_number <- function(value, name, minimum) {
  if (!is.numeric(value) || length(value) != 1 || is.na(value) ||
      !is.finite(value) || value != round(value) || value < minimum) {
    stop("Argument '", name, "' must be a single whole number of at least ",
         minimum, ".", call. = FALSE)
  }
  return(invisible(value))
}

# Whether the sampler reports its progress. NULL for silence, which is how a
# model that was not watched says so and keeps `model` free of a field that
# means nothing.
#
# The reporter behind this has been in the package since the C++ core arrived
# -- percentages, throttled to one line per percent -- with nothing to turn it
# on. A chain of 20000 draws over 196 series takes minutes and said nothing.
.check_verbose <- function(verbose) {
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("Argument 'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!verbose) {
    return(NULL)
  }
  return(TRUE)
}

# A time series as a named matrix, with its tsp intact.
#
# A univariate series arrives as a vector and everything downstream wants a
# matrix, so the conversion has to be assigned back: naming a vector is an
# error rather than a no-op, `as.matrix()` drops the `tsp`, hence the copy
# either side of it, and `scale()` would drop a vector's while keeping a
# matrix's. create_dfmodel() has always done this; create_favarmodel() did not,
# so a plain `ts` vector passed as its observed block reached the binding as a
# vector, came back as Rcpp's "Not a matrix." and was turned into `error = TRUE`
# -- after which the only thing the caller saw was that there were no draws.
#
# Keyed on the column names rather than on `dimnames` as a whole, so that a
# matrix which has row names but no column names is named here as well instead
# of carrying empty names through to the output. A `ts` matrix built without
# names already has R's own "Series 1" and is left alone.
.as_named_ts_matrix <- function(x, prefix) {

  if (!is.null(colnames(x))) {
    return(x)
  }

  tsp_temp <- stats::tsp(x)
  x <- stats::ts(as.matrix(x))
  stats::tsp(x) <- tsp_temp
  colnames(x) <- if (NCOL(x) == 1) prefix else paste0(prefix, seq_len(NCOL(x)))

  return(x)
}
