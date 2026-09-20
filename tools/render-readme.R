## Regenerates README.md from README.Rmd. Run from the package root:
##
##     Rscript tools/render-readme.R
##
## devtools::build_readme() is the usual route and produces the same thing, but
## it needs pandoc. This does not: the body of README.Rmd is already
## GitHub-flavoured markdown, so knitting it and dropping the YAML header knitr
## leaves behind is equivalent. The output_format in the Rmd stays
## github_document all the same, for whoever does have pandoc.

for (pkg in c("knitr", "dfmtools")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is needed to render the README.", call. = FALSE)
  }
}

if (!file.exists("README.Rmd")) {
  stop("Run this from the package root: README.Rmd is not here.", call. = FALSE)
}

## A chunk that fails must fail the render. knitr's default is to record the
## error in the output and carry on, which is how a README whose every forecast
## chunk was an error message once reached the front page of the repository.
## vignettes/precompile.R has set this for the same reason for as long as it has
## existed. Chunks that mean to show an error still set error = TRUE themselves.
knitr::opts_chunk$set(error = FALSE)

## R's own messages otherwise come out in the renderer's language, and a README
## is read in English whoever renders it. Set for this process only; the script
## is run with Rscript and takes nothing back to the session.
Sys.setenv(LANGUAGE = "en")

## The chunks call `library(dfmtools)`, so what is baked into README.md is the
## *installed* package, not the sources in this working tree. Install first:
##
##     R CMD INSTALL .
knitr::knit("README.Rmd", "README.md", encoding = "UTF-8")

lines <- readLines("README.md", encoding = "UTF-8", warn = FALSE)

if (length(lines) && lines[1] == "---") {
  closing <- which(lines == "---")
  if (length(closing) >= 2) {
    lines <- lines[-seq_len(closing[2])]
  }
}

while (length(lines) && !nzchar(lines[1])) {
  lines <- lines[-1]
}

## useBytes: readLines above marked these as UTF-8, and writeLines would
## otherwise re-encode them to the native codepage on Windows.
writeLines(lines, "README.md", useBytes = TRUE)

message("Wrote README.md (", length(lines), " lines).")
