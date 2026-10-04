# Deterministic terms of a factor model.
#
# A constant, a linear trend and seasonal dummies, built the way bvartools'
# create_bvarmodel() builds them for a VAR -- the same names, "const", "trend"
# and "season.1" onwards, the same trend starting at one in the first period,
# and the same dummies, "season.i" being one in the i-th period of the year and
# the last period of the year the one left out -- so that a model specified
# with the same arguments in either package carries the same columns.
#
# Where they enter is where the two packages part company. In a VAR they are
# regressors like any other. In a factor model they enter the measurement,
# x_t = lambda f_t + C d_t + u_t, so that the factors are deviations from them,
# and a FAVAR's observed variables deviate from C_obs d_t on top of that. See
# the BayesTS core's DfmNormalGammaInput and FavarNormalWishartInput.


# The two arguments as create_dfmodel() and create_favarmodel() take them, with
# bvartools' rule that seasonal dummies need a constant beside them: without one
# the dummies would be read as the levels of the seasons rather than as their
# differences from the one left out.
.check_deterministic_spec <- function(deterministic, seasonal) {
  if (!is.character(deterministic) || length(deterministic) != 1 ||
      !deterministic %in% c("none", "const", "trend", "both")) {
    stop("Argument 'deterministic' must be one of 'none', 'const', 'trend' or 'both'.",
         call. = FALSE)
  }
  if (!is.logical(seasonal) || length(seasonal) != 1 || is.na(seasonal)) {
    stop("Argument 'seasonal' must be TRUE or FALSE.", call. = FALSE)
  }
  if (seasonal && !deterministic %in% c("const", "both")) {
    stop("Argument 'deterministic' must be either 'const' or 'both' when using ",
         "'seasonal = TRUE'.", call. = FALSE)
  }
  invisible(TRUE)
}


# The deterministic terms over the periods of `x`, a time series object, as a
# time series matrix with one named column per term, or NULL for a model
# without any. `x` is only read for its time attributes.
.deterministic_terms <- function(x, deterministic, seasonal) {

  tt <- NROW(x)
  freq <- stats::frequency(x)
  terms <- NULL

  if (deterministic %in% c("const", "both")) {
    terms <- cbind(terms, const = rep(1, tt))
  }
  if (deterministic %in% c("trend", "both")) {
    terms <- cbind(terms, trend = seq_len(tt))
  }
  if (seasonal) {
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.",
              call. = FALSE)
    } else {
      period <- stats::cycle(stats::ts(seq_len(tt), start = stats::start(x), frequency = freq))
      for (i in seq_len(freq - 1)) {
        terms <- cbind(terms, as.numeric(period == i))
        colnames(terms)[ncol(terms)] <- paste0("season.", i)
      }
    }
  }

  if (is.null(terms)) {
    return(NULL)
  }

  stats::ts(terms, start = stats::start(x), frequency = freq)
}


# The deterministic terms of a model over the `n_ahead` periods after its
# sample, continued by the rule each was built by: a constant stays one, a trend
# goes on counting, and a seasonal dummy follows the calendar. Built by
# .deterministic_terms() over the sample and the horizon together, so that the
# horizon cannot disagree with the sample about which period is which.
.continue_deterministic <- function(object, n_ahead) {

  train <- object[["data"]][["deterministic"]]
  names <- colnames(train)
  known <- c("const", "trend", paste0("season.", seq_len(max(1, stats::frequency(train) - 1))))
  if (!all(names %in% known)) {
    stop("Could not identify all deterministic terms. Please specify argument ",
         "'deterministic' instead.", call. = FALSE)
  }

  deterministic <- if ("trend" %in% names) {
    if ("const" %in% names) "both" else "trend"
  } else if ("const" %in% names) "const" else "none"
  seasonal <- any(startsWith(names, "season."))

  tt <- NROW(train)
  span <- stats::ts(seq_len(tt + n_ahead), start = stats::start(train),
                    frequency = stats::frequency(train))
  terms <- .deterministic_terms(span, deterministic, seasonal)
  terms[tt + seq_len(n_ahead), names, drop = FALSE]
}


# The deterministic terms of the `n_ahead` forecast periods, h rows by one
# column per term: continued from the sample where `deterministic` is NULL, and
# otherwise taken from it, a time series that has to cover those periods, as in
# bvartools' prepare_forecast_input(). NULL for a model without any.
.forecast_deterministic <- function(object, n_ahead, deterministic = NULL) {

  train <- object[["data"]][["deterministic"]]
  if (is.null(train)) {
    if (!is.null(deterministic)) {
      stop("Argument 'deterministic' was given, but the model has no deterministic terms.",
           call. = FALSE)
    }
    return(NULL)
  }

  if (is.null(deterministic)) {
    return(unclass(.continue_deterministic(object, n_ahead)))
  }

  if (!"ts" %in% class(deterministic)) {
    stop("Argument 'deterministic' must be of class 'ts'.", call. = FALSE)
  }
  deterministic <- .as_named_ts_matrix(deterministic, "deterministic")
  names <- colnames(train)
  if (NCOL(deterministic) != length(names)) {
    stop("Argument 'deterministic' must have one column per deterministic term of the ",
         "model, ", length(names), " of them (", paste(names, collapse = ", "), ").",
         call. = FALSE)
  }

  freq <- stats::frequency(train)
  first <- stats::tsp(train)[2] + 1 / freq
  last <- first + (n_ahead - 1) / freq
  tsp_d <- stats::tsp(deterministic)
  tol <- 0.1 / freq
  if (abs(tsp_d[3] - freq) > 1e-8 || tsp_d[1] > first + tol || tsp_d[2] < last - tol) {
    stop("Argument 'deterministic' must cover the ", n_ahead, " forecast period(s), from ",
         .period_label(first, freq), " to ", .period_label(last, freq),
         ", at the frequency of the model (", freq, ").", call. = FALSE)
  }

  # By time rather than through stats::window(), which warns about a start or
  # an end it rounds to the nearest period.
  times <- as.numeric(stats::time(deterministic))
  rows <- vapply(first + (seq_len(n_ahead) - 1) / freq,
                 function(t) which(abs(times - t) < tol)[1], integer(1))
  x <- unclass(as.matrix(deterministic))[rows, , drop = FALSE]
  x <- matrix(as.numeric(x), n_ahead, length(names), dimnames = list(NULL, names))
  if (!all(is.finite(x))) {
    stop("Argument 'deterministic' must be finite over the forecast periods.", call. = FALSE)
  }
  x
}


# A period of a time series as stats::window() and stats::ts() take it, e.g.
# c(2019, 4), which is how a message asks for the periods a series has to cover.
.period_label <- function(time, frequency) {
  year <- floor(time + 1e-6)
  if (frequency == 1) {
    return(as.character(year))
  }
  paste0("c(", year, ", ", round((time - year) * frequency) + 1, ")")
}


# The normal prior on vec(C), `rows` by the number of terms, column by column:
# a zero mean and `spec$vinv` on the diagonal. Diagonal by construction, which
# is what the sampler needs -- it draws C a row at a time and refuses a
# precision coupling two rows.
.deterministic_prior <- function(spec, rows, n_det, name) {
  .check_coefficient_prior(spec, name)
  list(mu = matrix(0, rows * n_det),
       vinv = diag(spec[["vinv"]], rows * n_det))
}


# Where the chain starts C: the least squares coefficients of each column of
# `y` (one row per period) on the deterministic terms `d`, as vec(C) of the
# rows-by-terms matrix the sampler draws.
#
# Not a draw from the prior, which is what every other block starts from. A
# persistent factor can carry part of a level or a trend for a small price in
# its own transition, so the constant and the trend of C trade off against the
# level of the factors along a ridge the Gibbs sampler walks slowly. A chain
# started at a prior draw -- or at zero -- puts the level into the factors at
# its first iteration and gives it back over thousands more; one started here
# begins near where the posterior puts it.
.deterministic_initial <- function(y, d) {
  y <- unclass(as.matrix(y))
  d <- unclass(as.matrix(d))
  coefficients <- qr.coef(qr(d), y)
  coefficients[is.na(coefficients)] <- 0
  matrix(t(coefficients))
}


# Column names for the draws of vec(C): "<series>.<term>", series fastest,
# which is the order vec() of a series-by-terms matrix puts them in.
.deterministic_draw_names <- function(series, terms) {
  paste(rep(series, length(terms)), rep(terms, each = length(series)), sep = ".")
}
