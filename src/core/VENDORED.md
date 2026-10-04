# Vendored BayesTS core

`src/core/` and `inst/include/bayests/` are a copy of part of the core layer of
[BayesTS](https://github.com/franzmohr/BayesTS): the sampler declaration in
`include/bayests/` and the numerics in `src/core/`. Nothing else of that
project is here -- the core deliberately links neither HDF5 nor HighFive,
prints nothing and reads no files, which is what makes it embeddable in an R
package at all.

The copy is **BayesTS `065ed77`**, the tip of upstream's
`factor-deterministic-terms` branch, off `main` at `82ac49c` and merged into it
as `6fa1c75`, like every branch vendored before it, with its commits unchanged
-- so the files here are those of upstream `main`. The branch adds deterministic terms to all five factor models (`c2c49ea` --
`inputs.h`, `results.h`, `spec.h`, the five sampler headers, `inputs.cpp`,
`dfm_support.h`, `favar_support.h`, `factor_score.h` and the five models'
sources; `065ed77`, two comments): a constant, a trend or seasonal dummies in
the measurement as `C d_t`, and a FAVAR's observed variables deviating from
`C_obs d_t`. It is the second refresh to widen what a factor model can be asked
for, and it does so without moving a draw of a model without deterministic
terms: the new blocks consume no random numbers when `n` is zero, and upstream
compared the fingerprints of all 155 fixtures before against after, none moved.
Before that the copy was `696dc7d`, the tip of upstream's `structural-qvar`
branch, merged into `main` as `82ac49c`, after the 0.3.0 release
(tagged on `18a86c2`). Past `7d7c6ca` the files reached here moved for reasons
none of which changes a factor model's draws: the `bvs` flat-prior diagnostic
and the SSVS refusals (`fa6dd89`, `f2abd4e`, `a85abe2` -- `priors.h`,
`inputs.cpp`, `model_support.h`); the non-centred random walks of
`VarTvpStochvol`, `VarTvpGamma`, `VecTvpStochvol` and `VecTvpGamma` (`59c495f`,
`0a2a0c0`, `95bacdb`, `cea124b` -- `results.h`, `inputs.cpp` and three additions to the stochastic
volatility files); the constant VECs' cointegration prior (`b8d6c2c` --
`priors.h`, `inputs.cpp`); the smoother skipping an identity transition
(`4736936` -- `kalman_durbin_koopman_2002.cpp`, exact); and the guarded
forecast root (`62bc455` -- `forecast_states.h`, whose `covariance_root()`
`favar_normal_wishart.cpp` now calls in place of its own copy, which differed
only in flooring a negative eigenvalue with `abs()` rather than at zero, a
case no fixture reaches); a doc comment moved back onto
`validate_tvp_block()` (`ca1c322` -- `inputs.cpp`, comment only); the refusal
of a NaN or an infinity in any input, by name (`31735c8` -- `model_support.h`'s
`require_finite()`, called from every `validate()` in `inputs.cpp` and from
`scored_horizons()` in `predictive_score.h`); and Doxygen parameter comments
(`fd7ece5` -- `chan_jeliazkov_2009.cpp`, `stochvol_mixture.h`,
`dfm_support.h`, comments only); the non-centred parameterisation of the two
time-varying factor models' random walks (`1e21f94` -- `noncentred_support.h`,
which is new here, `dfm_tvp_gamma.cpp`, `dfm_tvp_stochvol.cpp`, `priors.h`,
`results.h`, `inputs.cpp`), which is the first refresh since this package
existed to widen what a factor model can be asked for rather than only what it
refuses; and `/model/n_iid`, an equation carrying no coefficients, which is a
VAR feature whose refusal every `validate()` calls (`3425586` -- `spec.h`,
`inputs.cpp`, `model_support.h`); and, past `d5ca82c`, panels not observed whole
and four options of a VAR's prior, and an i.i.d. block beside them (`567aa57`
to `2f9f010` -- `data.h`,
`inputs.h`, `priors.h`, `results.h`, `spec.h`, `spec.cpp`, `inputs.cpp`, and
new here `constraint_support.h`, `completion_support.h`, `shrinkage_support.h`,
`steady_state_support.h`, `constrained_var_path.{h,cpp}` and
`gig_hormann_leydold_2014.{h,cpp}`, which the
`validate()`s in `inputs.cpp` now reach); and the structural quantile VAR
(`2e4fe3e`, `696dc7d` -- `spec.h`, `results.h`, `inputs.cpp`), a grid of
quantiles on `VarNormalAld` whose numerics, `quantile_grid.h` and
`var_normal_ald.cpp`, nothing here reaches. Every factor model refuses the
constraint datasets and none reads the prior options, so nothing of it reaches
a factor model's draws; upstream verified the draws unchanged over 145
fixtures, then 148, 151, 153 and 154, the factor models' among them. The refreshes before sat at `696dc7d`, `2f9f010`, `afeb326`, `89b0495`, `638355a`, `d5ca82c`, `d92e581`,
`f8b42a1`, `94f81de`, `cea124b`, `ffbd208` and `4a64082`. 0.3.0 is not archived yet,
so there is no version DOI to name here; the concept DOI
<https://doi.org/10.5281/zenodo.22722531> resolves to the newest release
whenever one is cut. The last archived release is 0.2.0,
<https://doi.org/10.5281/zenodo.22765348>, which this copy is past.

Upstream commits that change nothing under `include/` or `src/core/` do not move
the copy off that commit; a refresh that copies anything newer has to update this
paragraph. Nothing enforces that -- the script compares the vendored *files* and
`inst/COPYRIGHTS`, never this prose -- and it has gone stale once already: the
refreshes that brought `forecast_states.h` and then the forecast score left the
paragraph claiming `v0.2.0` while the files were a release's worth of commits
past it. So check it against the upstream `git log` on the way out, as part of
the refresh rather than after it.

The same arrangement bvartools uses, and for the same reason: there is one
implementation of each sampler, upstream, and the R packages are translation
layers over it. `src/DfmNormalGamma.cpp`, `src/DfmNormalStochvol.cpp`,
`src/DfmTvpGamma.cpp`, `src/DfmTvpStochvol.cpp` and
`src/FavarNormalWishart.cpp` are the whole of this package's half of that --
each turns a `dfmodel` or `favarmodel` list into the corresponding
`bayests::...Input` and the draws back into a list. What the four dynamic factor
bindings share is in `src/dfm_r_translation.h`, which is one thing: the
permutation between R's `lower.tri()` ordering of the free loadings and the
row-major ordering the core draws them in. The FAVAR binding does not use it --
its identification is an identity block rather than a unit lower triangle, so
every row after the block is free across its whole width and the free elements
are a plain rectangle.

## Which files, and why not all of them

Only what the five factor models reach. This package has those five samplers;
upstream has twenty-two, and compiling the other seventeen would cost a minute of
build time and a larger shared object for code that is never called. bvartools
mirrors the whole core because it uses fifteen of them, which is a different
trade.

**`FavarNormalWishart` is now taken.** A factor augmented VAR is a factor model
-- it reaches `dfm_support.h` and the band sampler like the four above it -- so
bvartools skips it as it skips them, and this is where it goes. Adding it was two
lines in `entry` below, and the closure pulled in `favar_support.h` and
`core/algorithms/wishart.{h,cpp}` on its own, exactly as this section predicted
before it was done: five files, from twenty-seven to thirty-two. It arrives with
a binding and an R entry point -- `create_favarmodel()` and the four `favarmodel`
methods -- so it is not the sampler-nothing-can-call this section otherwise
declines.

Two things about it belong here rather than in the R docs. Its observed factors
go into `TrainData::f_obs`, not into `x` or `z`: those are the regressor layouts,
and the observed block is half of the *state*, appearing on the left of the
transition as well as the right. And its forecast is the one here that is wider
than `k` -- each horizon carries the panel followed by the observed factors,
which have no other object to go in.

Also here from that work is the second entry point on the band
sampler. `core/algorithms/chan_jeliazkov_2009.cpp` now carries
`chan_jeliazkov_2009_conditional`, which holds the trailing elements of every
state column at observed values instead of drawing them -- what a state vector
that is part data needs -- and `dfm_support.h` carries
`draw_conditional_factor_path()` beside `draw_factor_path()`, on the same
one-block-or-stack contract and with the same shift convention. Neither is
reached by the four samplers here. Both arrived through files this package
already vendors, and upstream verified the refactor that introduced them left
every fingerprint unmoved.

Each sampler after the first has been cheap to add, because they share
`dfm_support.h`, `model_support.h` and `chan_jeliazkov_2009`.
`DfmNormalStochvol` grew the closure by its own two files plus the stochastic
volatility mixture -- five in all, from sixteen files to twenty-one.
`DfmTvpGamma` grew it by four: its own two, and
`kalman_durbin_koopman_2002.{h,cpp}`, which it needs for the loading paths and
the transition path and which no constant-coefficient model reaches.
`DfmTvpStochvol` grew it by two, its own, having nothing left to reach for that
the three before it had not already brought. `FavarNormalWishart` grew it by
five: its own two, `favar_support.h`, and `wishart.{h,cpp}`, which no factor
model with a diagonal transition precision reaches. Thirty-two at that point.

Four have arrived since, none of them a sampler: `forecast_states.h` with the
refresh that carries a model's drift over the forecast horizon;
`factor_score.h` and `predictive_score.h` with `predictive_log_density()`, which
is how a forecast is scored against what its periods realised; and
`noncentred_support.h` with the refresh that put `DfmTvpGamma` and
`DfmTvpStochvol` on the non-centred parameterisation of their random walks.
**Thirty-six now**, and that is the number `inst/COPYRIGHTS` lists and the
refresh script counts.

The set is *computed* rather than listed. `tools/update-bayests-core.R` starts
from

    include/bayests/dfm_normal_gamma.h
    src/core/models/dfm_normal_gamma.cpp
    include/bayests/dfm_normal_stochvol.h
    src/core/models/dfm_normal_stochvol.cpp
    include/bayests/dfm_tvp_gamma.h
    src/core/models/dfm_tvp_gamma.cpp
    include/bayests/dfm_tvp_stochvol.h
    src/core/models/dfm_tvp_stochvol.cpp
    include/bayests/favar_normal_wishart.h
    src/core/models/favar_normal_wishart.cpp
    src/core/spec.cpp
    src/core/inputs.cpp

and follows every `#include` of a `"bayests/..."` or `"core/..."` header until
nothing new turns up, taking each header's `.cpp` alongside it. A hand-written
list would go stale the first time upstream added an include; this cannot, and
a file that stops being reachable is reported rather than left behind.

`spec.cpp` and `inputs.cpp` are named as entry points rather than discovered,
because nothing includes them -- they hold `VarSpec::validate()` and every
`Input::validate()`, which the samplers call across translation units. Each
sampler has to be named for the same reason as the others: none is reachable from
another's includes, so a sampler left off that list is simply not vendored, and
the omission shows up as a link error rather than as a warning.

| Upstream | Here |
| --- | --- |
| `include/bayests/*.h` | `inst/include/bayests/*.h` |
| `src/core/{spec,inputs}.cpp` | `src/core/` |
| `src/core/models/*` | `src/core/models/` |
| `src/core/algorithms/*` | `src/core/algorithms/` |

A refresh is

```bash
Rscript tools/update-bayests-core.R path/to/BayesTS/source_tree
```

It reports what it changed, what it skipped and anything here that upstream no
longer reaches, and it is idempotent -- running it twice changes nothing the
second time.

The copied files keep their `SPDX-License-Identifier: BSD-3-Clause` headers.
That is compatible with this package's GPL (>= 2) and the copyright holder is
the same person, so there is nothing to reconcile. `inst/COPYRIGHTS` spells the
difference out for a CRAN reviewer and lists the files by name. That list has to
match what is actually vendored, and the refresh script compares the two and
refuses to finish when they disagree -- it does not edit the list, so a refresh
that added or removed a file stops until the list is brought into line by hand.

## Local modifications

Keep this list short; every entry is something a refresh has to reapply.

1. **Armadillo comes from RcppArmadillo.** Upstream opens with

   ```cpp
   #define ARMA_DONT_PRINT_FAST_MATH_WARNING
   #include <armadillo>
   ```

   and every such `#include <armadillo>` is replaced by
   `#include "bayests/arma.h"`. That header is upstream's own and resolves to
   whatever `BAYESTS_ARMA_HEADER` names; `src/Makevars` defines it as
   `<RcppArmadillo.h>`, which is what points Armadillo's RNG at R's, so
   `set.seed()` reaches the sampler. The rule is uniform -- every copied file
   that includes Armadillo, header or source -- which is why it can be applied
   by a script rather than kept as a list of file names that goes stale.

   Losing this patch is the only way to break the vendored core silently: it
   compiles, links and runs, and merely stops honouring `set.seed()`. Two things
   guard it. The refresh script refuses to finish if any copied file still
   reaches for `<armadillo>`, and `src/bayests_rng_guard.cpp` turns the
   remaining case -- a file that never had the include but ends up seeing the
   wrong Armadillo anyway -- into a compile error, by asserting that
   `arma::arma_rng::rng_method` came out as R's.

2. **Nothing else.** No source is edited on the way in. A change the sampler
   needs is a change to make upstream.

## What is not vendored

`inst/include/bayests_reporter.h`, `src/bayests_r_io.h` and
`src/dfm_r_translation.h` belong to this package. The first implements the core's
`Reporter` contract against R's console and `Rcpp::checkUserInterrupt()`; the
second translates between an `Rcpp::List` and the core's structs; the third holds
the one piece of that translation the four dynamic factor bindings need and the
FAVAR binding does not. All
three are GPL (>= 2) like the rest of the package, and `inst/COPYRIGHTS` says
so.
