#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.update_Ha)]]
arma::sp_mat update_Ha(arma::mat& a, int& n_factors, int& p, int& tt) {

  // Transform matrix
  arma::mat amat = arma::reshape(a, n_factors, n_factors * p);
  arma::mat amat_t = arma::zeros<arma::mat>(n_factors * p, n_factors);
  for (int i = 0; i < p; i++) {
    amat_t.rows(i * n_factors, (i + 1) * n_factors - 1) = -amat.cols(i * n_factors, (i + 1) * n_factors - 1);
  }

  // Update H_a
  arma::sp_mat result = arma::speye<arma::sp_mat>(tt * n_factors, tt * n_factors);

  for (int i = 0; i < (tt - 1); i++) {
    if (i < tt - p) {
      result.submat((i + 1) * n_factors, i * n_factors,
                    (i + 1) * n_factors + n_factors * p - 1, (i + 1) * n_factors - 1) = amat_t;
    } else {
      result.submat((i + 1) * n_factors, i * n_factors,
                    (i + 1) * n_factors + n_factors * (p - (tt - i - 1)) - 1, (i + 1) * n_factors - 1) = amat_t.rows(0, n_factors * (p - (tt - i - 1)) - 1);
    }
  }

  return result;
}
