# CRAN submission status

**Not ready to submit yet.** One thing blocks it, and a second change follows
from it; this file tracks both so the submission itself is routine once they're
done.

1. `bvartools (>= 0.3.0.9000)` is a hard `Depends`, and that is a development
   version: CRAN has 0.3.0 (published 2026-09-11). CRAN resolves dependencies
   from CRAN itself, so submitting dfmtools before a bvartools release at or
   above that version is accepted will fail outright. Submit bvartools first,
   and raise the requirement here to whatever it is released as.
2. `DESCRIPTION` carries `Remotes: franzmohr/bvartools` so that
   `remotes::install_github()` and the GitHub Actions check pick bvartools up
   from GitHub in the meantime. Drop that field before the actual submission
   -- once the required bvartools is on CRAN it is no longer needed, and CRAN
   flags it as an unknown field otherwise.

The URL check no longer flags anything: `https://github.com/franzmohr/dfmtools`
is public, so the `URL`, `BugReports` and `inst/CITATION` links resolve, and the
Buy Me a Coffee badge in the README links to `https://buymeacoffee.com/franzmohr`
rather than the `www.` address that redirects there.

## Test environments

* Local: Windows 11 x64 (build 26200), R 4.6.1 (ucrt), Rtools45
  (g++ 14.2.0, x86_64-w64-mingw32), `R CMD check --as-cran` on the tarball
  `R CMD build` makes from the `v0.1.0` tag.
* Not yet run: win-builder, R-hub, macOS. Do at least one of these before
  submitting -- the vendored C++17 core (`src/core/`) has only been built
  and exercised on Windows/MinGW so far.

`Depends` says `R (>= 4.0.0)` rather than the `(>= 3.5)` it used to, and rather
than the `(>= 3.5)` bvartools says. `src/Makevars` sets `CXX_STD = CXX17`, and
3.5 predates a toolchain that reliably honours it; the oldest R the check
matrix exercises is `oldrel-1`, so nothing ever tested the old claim. Stricter
and true beats looser and untested, but it does mean this package and its own
dependency now disagree about their floor -- worth reverting in one line if
they should stay in step.

## R CMD check results

0 errors | 0 warnings | 3 notes

* `checking installed package size ... NOTE`: the installed package is about
  24Mb, nearly all of it `libs`. That is the vendored BayesTS core, which is
  five Gibbs samplers of templated C++ compiled with debugging symbols by
  default; it is the size the package is and there is nothing to strip that
  CRAN does not strip itself. Recent R reports this as `INFO` rather than
  `NOTE`, so a local check may not show it where CRAN's does.
* `checking CRAN incoming feasibility ... NOTE`, which carries two items:
    * `New submission` -- expected, this is a first submission.
    * `Unknown, possibly misspelled, fields in DESCRIPTION: 'Remotes'` -- goes
      away once the `Remotes:` field above is dropped.
* `checking HTML version of manual ... NOTE`: `Skipping checking math
  rendering: package 'V8' unavailable` -- V8 isn't installed in this local
  environment. An artifact of this machine, not the package; it does not
  appear on CRAN's own check machines.

Pandoc isn't on this machine's `PATH` either. The check above put the copy
bundled with Positron (`resources/app/quarto/bin/tools/pandoc.exe`) on it; a
check run without that adds a third local-only note, `Files 'README.md' or
'NEWS.md' cannot be checked without 'pandoc' being installed`.

The tests pass: `[ FAIL 0 | WARN 0 | SKIP 4 | PASS 910 ]`. The four skips are
the three cases of `test-parallel.R` that start worker processes, and the unload
test of `test-package-hooks.R`, which starts one of its own. All four run here:
with `NOT_CRAN` set the suite skips nothing against the installed package, which
is worth doing before a submission because the skipped four are the ones CRAN
never exercises.

## Downstream dependencies

None -- this is a new package.
