# The package's own load and unload hooks.
#
# .onUnload() is the one piece of R in this package that nothing else calls,
# so nothing else notices when it is wrong -- and what it would be wrong about
# is the name of the shared object, which only shows up as a warning at the end
# of a session or a failed devtools::load_all() cycle.

test_that("the shared object is registered under the name .onUnload uses", {

  loaded <- getLoadedDLLs()
  expect_true("dfmtools" %in% names(loaded))

  # library.dynam.unload("dfmtools", libpath) finds the object by this name, so
  # the two have to agree. .registration = TRUE is what puts the entry points
  # in the same place.
  expect_true(is.loaded("_dfmtools_DfmNormalGammaCoefficients"))
  expect_true(is.loaded("_dfmtools_FavarNormalWishartCoefficients"))
})

test_that("the namespace unloads cleanly", {

  # In a process of its own: unloading the package under test would take the
  # rest of the suite with it. The check itself does this from R CMD check's
  # "can be unloaded cleanly" step, but only for the namespace as installed --
  # this covers .onUnload() running as part of it.
  skip_on_cran()
  skip_on_os("solaris")

  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  skip_if_not(file.exists(rscript), "no Rscript to start a second process with")

  code <- paste(
    'library(dfmtools)',
    'stopifnot("dfmtools" %in% names(getLoadedDLLs()))',
    'unloadNamespace("dfmtools")',
    'stopifnot(!("dfmtools" %in% names(getLoadedDLLs())))',
    'cat("unloaded\n")',
    sep = "; ")

  out <- suppressWarnings(system2(rscript,
                                  c("--vanilla", "-e", shQuote(code)),
                                  stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("unloaded", out, fixed = TRUE)),
              label = paste(out, collapse = " | "))
})
