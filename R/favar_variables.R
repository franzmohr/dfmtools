#' Resolve a variable to its position
#'
#' Shared by the loading restrictions of \code{\link{add_priors.favarmodel}} and
#' the lookups of \code{\link{irf.favarmodel}}, so that a name means the same
#' thing and a wrong one reads the same way in both.
#'
#' @param value the name or index provided by the user.
#' @param names the names it is looked up among.
#' @param argument the argument's name, for the error message.
#'
#' @return An integer position.
#'
#' @noRd
.favar_position <- function(value, names, argument) {

  if (is.null(value)) {
    stop("Argument '", argument, "' must be provided. Available are: ",
         paste(names, collapse = ", "), ".")
  }
  if (length(value) != 1) {
    stop("Argument '", argument, "' must be a single variable.")
  }

  if (is.character(value)) {
    position <- which(names == value)
    if (length(position) == 0) {
      stop("Variable '", value, "' is not available as '", argument,
           "'. Available are: ", paste(names, collapse = ", "), ".")
    }
    # A panel column and an observed one may carry the same name, and which of
    # the two was meant is not something to guess at.
    if (length(position) > 1) {
      stop("Name '", value, "' is not unique among the variables available as '",
           argument, "'. Use its index, ", paste(position, collapse = " or "),
           ", instead.")
    }
    return(position)
  }

  position <- suppressWarnings(as.integer(value))
  if (is.na(position) || position < 1 || position > length(names)) {
    stop("Argument '", argument, "' must be between 1 and ", length(names), ".")
  }

  return(position)
}

#' The panel's column names
#'
#' @param object an object of class \code{'favarmodel'}.
#' @param m the number of panel series.
#'
#' @return A character vector of length \code{m}.
#'
#' @noRd
.favar_panel_names <- function(object, m) {

  names <- colnames(object[["data"]][["x"]])
  if (is.null(names)) {
    names <- paste0("x", seq_len(m))
  }

  return(names)
}
