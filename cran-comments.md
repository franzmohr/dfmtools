# CRAN submission status

**Not ready to submit yet.** One thing blocks it, and a second change follows
from it; this file tracks both so the submission itself is routine once they're
done.

1. `bvartools (>= 1.0.0)` is a hard `Depends`, and CRAN has 0.3.0 (published
   2026-09-11). CRAN resolves dependencies from CRAN itself, so submitting
   dfmtools before bvartools 1.0.0 is accepted will fail outright. Submit
   bvartools first.
2. `DESCRIPTION` carries `Remotes: franzmohr/bvartools` so that
   `remotes::install_github()` and the GitHub Actions check pick bvartools up
   from GitHub in the meantime. Drop that field before the actual submission
   -- once bvartools 1.0.0 is on CRAN it is no longer needed, and CRAN flags
   it as an unknown field otherwise.

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

## R CMD check results

0 errors | 0 warnings | 2 notes

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

The tests pass: `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 715 ]`. The one skip is
`test-agent-docs.R`, which skips itself on CRAN.

## Downstream dependencies

None -- this is a new package.
