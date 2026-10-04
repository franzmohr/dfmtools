#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "dfm_r_translation.h"

// The permutation between R's ordering of a dynamic factor model's free
// loadings and the core's, for the R code that writes a model file.
//
// A model file is in the core's ordering, since the BayesTS command line reads
// it into the same structs these bindings fill. Exporting the bindings' own
// permutation rather than writing it again in R keeps the two ways into the
// core -- through memory and through a file -- from disagreeing: a permuted
// prior is still a prior, and nothing would notice.
//
// One-based, as R indexes: element `i` is the R position of the core's `i`-th
// free loading, so `core <- r[order]` and `r[order] <- core`.

// [[Rcpp::export(.lambda_row_major_order)]]
Rcpp::IntegerVector lambda_row_major_order_r(const int m, const int n) {

  const arma::uvec order = dfmtools::lambda_row_major_order(m, n);

  Rcpp::IntegerVector result(order.n_elem);
  for (arma::uword i = 0; i < order.n_elem; i++) {
    result[i] = static_cast<int>(order(i)) + 1;
  }
  return result;
}
