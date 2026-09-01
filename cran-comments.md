# CRAN submission status

**Not ready to submit yet.** Two things need to happen first; this file
tracks them so the submission itself is routine once they're done.

1. `bvartools (>= 1.0.0)` is a hard `Depends`, and only 0.2.4 is on CRAN.
   CRAN resolves dependencies from CRAN itself, so submitting dfmtools before
   bvartools 1.0.0 is accepted will fail outright. Submit bvartools first.
2. `DESCRIPTION` carries `Remotes: franzmohr/bvartools` so that
   `remotes::install_github()` and the GitHub Actions check pick bvartools up
   from GitHub in the meantime. Drop that field before the actual submission
   -- once bvartools 1.0.0 is on CRAN it is no longer needed, and CRAN flags
   it as an unknown field otherwise.
3. `DESCRIPTION`'s `URL` and `BugReports`, and the same links in
   `man/dfmtools-package.Rd` and `inst/CITATION`, point at
   `https://github.com/franzmohr/dfmtools`, which currently 404s -- the
   local checkout has no `git remote` configured. Push the repository (or
   make it public) before submitting, or CRAN's URL check will flag it.

## Test environments

* Local: Windows 11, R 4.6.1, Rtools45 (x86_64-w64-mingw32), `R CMD check --as-cran`
* Not yet run: win-builder, R-hub, macOS. Do at least one of these before
  submitting -- the vendored C++17 core (`src/core/`) has only been built
  and exercised on Windows/MinGW so far.

## R CMD check results

0 errors | 0 warnings | 3 notes locally, of which one is expected to persist:

* `Unknown, possibly misspelled, fields in DESCRIPTION: 'Remotes'` -- goes
  away once the `Remotes:` field above is dropped.
* `Found the following (possibly) invalid URLs` (both `github.com/franzmohr/dfmtools`
  links, 404) -- goes away once the repository is pushed/public.
* `New submission` -- expected, this is a first submission.

Two further notes are artifacts of this machine, not the package, and do not
appear on CRAN's own check machines:

* `Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
  installed` -- pandoc isn't installed in this local environment.
* `Skipping checking math rendering: package 'V8' unavailable` -- V8 isn't
  installed in this local environment.

## Downstream dependencies

None -- this is a new package.
