// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef DFMTOOLS_DFM_R_TRANSLATION_H
#define DFMTOOLS_DFM_R_TRANSLATION_H

#include <RcppArmadillo.h>

#include <algorithm>

// What both dynamic factor model bindings have to translate, as opposed to what
// either of them reads on its own.
//
// There is one thing in that category and it is the ordering of the free
// loadings. It sat in DfmNormalGamma.cpp until DfmNormalStochvol.cpp wanted it
// too; a second copy of a permutation is the kind of thing that goes wrong in
// one place only and stays wrong, since a permuted prior is a prior.
namespace dfmtools
{

/// R keeps the free loadings in `lower.tri()` order -- column-major, so every
/// free row of column 0, then of column 1, and so on -- while the core wants
/// them row by row, which is the order its equation-by-equation draw consumes
/// them in and what makes each equation's slice of the prior contiguous.
///
/// This is that permutation: element `i` of the result is the R position of the
/// core's `i`-th free loading, so `core = r.elem(order)` and, for the prior
/// precision, `core = r.submat(order, order)`.
///
/// It matters for the prior, not for the starting value alone. `add_priors()`
/// builds that precision as `diag(vinv, n_lambda)`, which reads the same either
/// way, so a caller who takes the default cannot tell; a caller who supplies
/// anything else can.
inline arma::uvec lambda_row_major_order(const int m, const int n)
{

  arma::umat r_pos(m, n, arma::fill::zeros);
  arma::uword pos = 0;
  for (int j = 0; j < n; j++) {
    for (int i = j + 1; i < m; i++) {
      r_pos(i, j) = pos++;
    }
  }

  arma::uvec order(pos);
  arma::uword core = 0;
  for (int i = 1; i < m; i++) {
    for (int j = 0; j < std::min(i, n); j++) {
      order(core++) = r_pos(i, j);
    }
  }
  return order;
}

} // namespace dfmtools

#endif // DFMTOOLS_DFM_R_TRANSLATION_H
