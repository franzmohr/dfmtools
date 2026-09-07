# Vendored BayesTS core

`src/core/` and `inst/include/bayests/` are a copy of part of the core layer of
[BayesTS](https://github.com/franzmohr/BayesTS): the sampler declaration in
`include/bayests/` and the numerics in `src/core/`. Nothing else of that
project is here -- the core deliberately links neither HDF5 nor HighFive,
prints nothing and reads no files, which is what makes it embeddable in an R
package at all.

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
upstream has twenty, and compiling the other fifteen would cost a minute of build
time and a larger shared object for code that is never called. bvartools mirrors
the whole core because it uses fourteen of them, which is a different trade.

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
model with a diagonal transition precision reaches. Thirty-two now.

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
