# Vendored BayesTS core

`src/core/` and `inst/include/bayests/` are a copy of part of the core layer of
[BayesTS](https://github.com/franzmohr/BayesTS): the sampler declaration in
`include/bayests/` and the numerics in `src/core/`. Nothing else of that
project is here -- the core deliberately links neither HDF5 nor HighFive,
prints nothing and reads no files, which is what makes it embeddable in an R
package at all.

The same arrangement bvartools uses, and for the same reason: there is one
implementation of each sampler, upstream, and the R packages are translation
layers over it. `src/DfmNormalGamma.cpp` and `src/DfmNormalStochvol.cpp` are the
whole of this package's half of that -- each turns a `dfmodel` list into the
corresponding `bayests::...Input` and the draws back into a list. What the two
share is in `src/dfm_r_translation.h`, which is one thing: the permutation
between R's `lower.tri()` ordering of the free loadings and the row-major
ordering the core draws them in.

## Which files, and why not all of them

Only what the two dynamic factor models reach. This package has those two
samplers; upstream has fifteen, and compiling the other thirteen would cost a
minute of build time and a larger shared object for code that is never called.
bvartools mirrors the whole core because it uses twelve of them, which is a
different trade.

The second sampler was cheap to add: it shares `dfm_support.h`,
`model_support.h` and `chan_jeliazkov_2009` with the first, so the closure grew
by its own two files plus the stochastic volatility mixture -- five in all, from
sixteen files to twenty-one.

The set is *computed* rather than listed. `tools/update-bayests-core.R` starts
from

    include/bayests/dfm_normal_gamma.h
    src/core/models/dfm_normal_gamma.cpp
    include/bayests/dfm_normal_stochvol.h
    src/core/models/dfm_normal_stochvol.cpp
    src/core/spec.cpp
    src/core/inputs.cpp

and follows every `#include` of a `"bayests/..."` or `"core/..."` header until
nothing new turns up, taking each header's `.cpp` alongside it. A hand-written
list would go stale the first time upstream added an include; this cannot, and
a file that stops being reachable is reported rather than left behind.

`spec.cpp` and `inputs.cpp` are named as entry points rather than discovered,
because nothing includes them -- they hold `VarSpec::validate()` and every
`Input::validate()`, which the samplers call across translation units. Each
sampler has to be named for the same reason as the other: neither is reachable
from the other's includes, so a sampler left off that list is simply not
vendored, and the omission shows up as a link error rather than as a warning.

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
the one piece of that translation both dynamic factor bindings need. All three
are GPL (>= 2) like the rest of the package, and `inst/COPYRIGHTS` says so.
