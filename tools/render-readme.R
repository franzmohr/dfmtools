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
