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

`bvartools (>= 0.3.0.9000)` is a development version and is not on CRAN. It
comes from GitHub:

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

* **The README and the vignettes are generated, and from the *installed*
  package.** Both run real samplers, so both bake in whatever `library(dfmtools)`
  loads rather than the working tree. `R CMD INSTALL .` first, then

  ```sh
  Rscript tools/render-readme.R                          # README.md
  Rscript vignettes/precompile.R                         # both vignettes
  Rscript vignettes/precompile.R favar-monetary-policy   # or just one
  ```

  Regenerate them whenever a change alters what they show -- an API that moved,
  an argument that was added, a draw that comes out differently. Skipping this
  is how a `README.md` full of error messages and a vignette teaching an API
  that no longer exists both reached the repository at once. Both renderers set
  `error = FALSE`, so a broken chunk now stops them instead of being written
  into the output.
* **Style.** Follow the surrounding code rather than a style guide: two-space
  indent, `<-` for assignment, `stats::`/`coda::` prefixes on imported
  functions, and comments that say why rather than what.
* **NEWS.md.** A user-visible change gets an entry under the development
  version.

## Releasing, and the DOI

Every release is archived on [Zenodo](https://zenodo.org), which mints a DOI for
it. Two kinds of DOI come out of that: a **version DOI** for the single release,
and a **concept DOI** that always resolves to the newest one. The concept DOI is
the one to cite and the one the badge points at.

**One-time setup.** On Zenodo, log in with the GitHub account, open *GitHub* in
the account menu and flip the switch next to `franzmohr/dfmtools`. This has to
happen **before** the release is created -- Zenodo only receives releases
published after the switch is on, and an earlier one has to be uploaded by hand.
The repository must be public for the switch to appear.

**Metadata.** `.zenodo.json` is what Zenodo reads for the deposit: title,
description, author and ORCID, license, keywords. It deliberately carries no
`version` field, so Zenodo takes the version from the release tag and there is
one fewer file to keep in step. Zenodo archives the repository as GitHub
delivers it, so `.Rbuildignore` has no bearing on what the archive contains --
`^\.zenodo\.json$` is there only to keep the file out of the R tarball.

**Cutting a release.**

1. Bump `Version` and `Date` in `DESCRIPTION`, and `version` and
   `date-released` in `CITATION.cff`, to the same values.
2. Move the `NEWS.md` entries from *development version* to the release version.
3. `devtools::check()` clean, and `README.md` and the vignettes regenerated
   from their sources if any of what they show changed -- see *Making a change*
   above for the two commands.
4. Tag and push: `git tag -a v0.1.0 -m "dfmtools 0.1.0" && git push origin v0.1.0`.
5. Publish a GitHub release for the tag -- publishing it, not just pushing the
   tag, is what notifies Zenodo.
6. Check the new deposit on Zenodo before relying on it, in particular that the
   license came through as GPL-2.0-or-later; the license vocabulary is the one
   field of `.zenodo.json` that fails quietly.

**After the first release only.** Zenodo now has a concept DOI. Put it in three
places, all of which have a comment marking the spot: `CITATION.cff`
(`identifiers:`), `inst/CITATION` (a `doi =` argument on the `bibentry()` and in
`textVersion`), and the badge row at the top of `README.Rmd` --

```
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

-- after which `README.md` is regenerated. Later releases need none of this
again: the concept DOI does not change.

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
