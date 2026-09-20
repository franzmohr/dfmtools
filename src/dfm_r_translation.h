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

/// A value put on the core's ordering where it is the size the permutation is
/// for, left alone where it is not, and left empty where there is none.
///
/// The middle case is the point. These used to be `if (size matches) assign`,
/// so a wrong-sized argument was quietly not passed on and the core's
/// validator -- which does reject it, by size -- reported the size of what it
/// had rather than of what the caller wrote: "prior precision of lambda must
/// be 7x7, got 0x0" for a 3 x 3. Handing the value over unpermuted keeps the
/// rejection and makes the message name the real size. Permuting it would not
/// be possible and would not help.
inline arma::vec permute_free_loadings(const arma::vec &value, const arma::uvec &order)
{
  if (value.n_elem == order.n_elem) {
    return arma::vec(value.elem(order));
  }
  return value;
}

inline arma::mat permute_free_loadings_rows(const arma::mat &value, const arma::uvec &order)
{
  if (value.n_rows == order.n_elem) {
    return arma::mat(value.rows(order));
  }
  return value;
}

inline arma::mat permute_free_loadings_both(const arma::mat &value, const arma::uvec &order)
{
  if (value.n_rows == order.n_elem && value.n_cols == order.n_elem) {
    return arma::mat(value.submat(order, order));
  }
  return value;
}

} // namespace dfmtools

#endif // DFMTOOLS_DFM_R_TRANSLATION_H
