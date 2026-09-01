# Contributing to dfmtools

Bug reports, questions and pull requests are all welcome. A couple of things
about this package are worth knowing before you start, because they decide
*where* a change belongs.

## Where the sampler lives

`src/core/` and `inst/include/bayests/` are **vendored copies** of part of the
core layer of [BayesTS](https://github.com/franzmohr/BayesTS). They are not
maintained here. A change to the numerics of the sampler is a change to make
upstream; a refresh then brings it in:

```bash
Rscript tools/update-bayests-core.R path/to/BayesTS/source_tree
```

`src/core/VENDORED.md` explains which files are copied, how that set is computed
from the `#include` graph, and the one modification applied on the way in
(Armadillo has to arrive as RcppArmadillo, or `set.seed()` silently stops
reaching the sampler). A pull request that edits a vendored file directly will
be undone by the next refresh, so it will be asked for upstream instead.

What *is* maintained here is the translation layer: `src/DfmNormalGamma.cpp`,
`src/bayests_r_io.h`, `inst/include/bayests_reporter.h` and everything in `R/`.

## Setting up

`bvartools (>= 1.0.0)` is not on CRAN yet. It comes from GitHub:

```r
# install.packages("remotes")
remotes::install_github("franzmohr/bvartools")
remotes::install_deps("path/to/dfmtools", dependencies = TRUE)
```

Building the package needs a C++17 toolchain and GNU make — Rtools on Windows,
Xcode command line tools on macOS.

## Making a change

* **Documentation is generated.** `man/*.Rd` and `NAMESPACE` come from the
  roxygen blocks in `R/`. Edit the `R/` file and run `devtools::document()`;
  never edit an `.Rd` by hand.
* **`src/RcppExports.*` are generated too**, by `Rcpp::compileAttributes()`.
* **Tests.** `tests/testthat/`, testthat edition 3. New behaviour comes with a
  test; a bug fix comes with a test that failed before it.
  `devtools::test()` runs them in a few seconds — the fixtures in
  `helper-dfm.R` simulate a 40-quarter, 4-variable model rather than use
  `bem_dfmdata`, which keeps the sampler fast enough to test properly.
* **Check before opening a pull request.** `R CMD check` should be clean:

  ```r
  devtools::check()
  ```

* **Style.** Follow the surrounding code rather than a style guide: two-space
  indent, `<-` for assignment, `stats::`/`coda::` prefixes on imported
  functions, and comments that say why rather than what.
* **NEWS.md.** A user-visible change gets an entry under the development
  version.

## Reporting a bug

Include a reproducible example — ideally built from `bem_dfmdata` or simulated
data in the report itself — plus the output of `sessionInfo()`. If the report is
about draws being wrong, say which parameter block and at what `n` and `p`; the
identifying restrictions on the loading matrix mean the interesting failures
tend to appear only above one factor.

## Licensing

Contributions are accepted under the package's license, GPL (>= 2). Note that
the vendored files carry BSD-3-Clause SPDX headers; `inst/COPYRIGHTS` lists them
and explains why both licenses are present.
