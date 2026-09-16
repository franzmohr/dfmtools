#' Seed of the Posterior Simulation
#'
#' Sets the seed with which the posterior of a dynamic factor model or a factor
#' augmented VAR is simulated.
#'
#' @param object an object of class \code{'dfmodel'} or \code{'favarmodel'},
#' usually after \code{\link{add_initial_values.dfmodel}} or
#' \code{\link{add_initial_values.favarmodel}}.
#' @param seed a non-negative whole number no larger than
#' \code{.Machine$integer.max}.
#' @param ... not used.
#'
#' @details
#' Calling this function is optional. \code{add_initial_values()} already stores
#' a seed, drawn from R's random number generator, in element \code{seed} of
#' \code{object$model}. \code{add_seed()} replaces it, for example to give a
#' model a seed that depends neither on the state of R's generator nor on the
#' worker of a cluster that happens to simulate it.
#'
#' \code{add_posterior_coefficients()} draws with that seed: R's generator is set
#' to it, with R's default kinds, for the simulation and put back as it was
#' afterwards. A model with a given seed therefore gives the same draws however
#' R's generator stands, and a call of \code{set.seed()} between
#' \code{add_initial_values()} and \code{add_posterior_coefficients()} does not
#' change them. \code{add_posterior_forecasts()} is not seeded this way and
#' draws from R's generator as it stands.
#'
#' The generic is that of \pkg{bvartools}, together with its methods for lists
#' of models; see \code{\link[bvartools]{add_seed}}. A 'modellist' returned by
#' \code{\link{create_dfmodel}} or \code{\link{create_favarmodel}} gets the
#' seeds \code{seed}, \code{seed + 1}, \ldots, one per model, as in the example.
#' \code{add_initial_values()} on such a list already gives each of its models a
#' seed of its own. Seeding the models of a list needs a \pkg{bvartools} whose
#' list method counts the models of other packages; earlier development versions
#' left them without a seed.
#'
#' @return \code{object} with the seed set.
#'
#' @examples
#' data("bem_dfmdata")
#'
#' model <- create_dfmodel(x = bem_dfmdata, p = 1, n = 1,
#'                         iterations = 20, burnin = 10)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#'
#' # add_initial_values() has set a seed; replace it
#' model <- add_seed(model, 20260916)
#' model$model$seed
#'
#' # The models of a list get the seeds 100 and 101
#' models <- create_dfmodel(x = bem_dfmdata, p = 1:2, n = 1,
#'                          iterations = 20, burnin = 10)
#' models <- add_priors(models)
#' models <- add_seed(models, 100)
#'
#' @export
add_seed.dfmodel <- function(object, seed, ...) {
  object[["model"]][["seed"]] <- .check_seed(seed)
  object
}

#' @rdname add_seed.dfmodel
#' @export
add_seed.favarmodel <- function(object, seed, ...) {
  object[["model"]][["seed"]] <- .check_seed(seed)
  object
}

# The three helpers below are those of bvartools' R/add_seed.R, which does not
# export them. They are copied rather than reached with ':::' so that a model
# here is seeded exactly as a VAR there, and should change only with them.

# A seed is stored as an integer. A plain 20260916 in R is a double.
.check_seed <- function(seed) {
  if (length(seed) != 1 || !is.numeric(seed) || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) {
    stop("Argument 'seed' must be a single whole number between 0 and ",
         .Machine$integer.max, ".", call. = FALSE)
  }
  as.integer(seed)
}

# The seed add_initial_values() gives a model that has none, from R's generator
# as it stands.
.draw_model_seed <- function() {
  sample.int(.Machine$integer.max, 1L)
}

# Evaluates 'expr' with R's generator set to 'seed', with R's default kinds so
# that a model draws the same on a cluster worker running L'Ecuyer-CMRG as in a
# plain session, and puts the generator back afterwards, kinds included. Without
# a seed 'expr' is evaluated with the generator as it stands. 'expr' is a
# promise and is only forced after set.seed().
#
# The samplers reach R's generator because the vendored core draws through
# RcppArmadillo; see src/core/VENDORED.md. That wiring is what makes this work at
# all.
.with_model_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  seed <- .check_seed(seed)

  global <- globalenv()
  # Before RNGkind(), which initialises the generator when it has no state yet.
  had_state <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_state) {
    old_state <- get(".Random.seed", envir = global, inherits = FALSE)
  }
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (had_state) {
      assign(".Random.seed", old_state, envir = global)
    } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  }, add = TRUE)

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  expr
}
