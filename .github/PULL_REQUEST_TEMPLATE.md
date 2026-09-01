<!--
Thanks for the pull request. .github/CONTRIBUTING.md has the details; the
checklist below is the short version.
-->

## What this changes

<!-- One or two sentences, and the issue number if there is one. -->

## Checklist

- [ ] `devtools::check()` is clean (no errors, warnings or new notes)
- [ ] `devtools::test()` passes, and new behaviour or a fixed bug comes with a test
- [ ] `devtools::document()` was run if any roxygen block changed — `man/` and
      `NAMESPACE` are generated, not edited
- [ ] `Rcpp::compileAttributes()` was run if any `// [[Rcpp::export]]` changed
- [ ] `NEWS.md` has an entry if the change is user-visible
- [ ] No file under `src/core/` or `inst/include/bayests/` is edited — those are
      vendored from [BayesTS](https://github.com/franzmohr/BayesTS) and a change
      there has to be made upstream (see `src/core/VENDORED.md`)
