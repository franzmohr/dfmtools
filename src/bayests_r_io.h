// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef DFMTOOLS_BAYESTS_R_IO_H
#define DFMTOOLS_BAYESTS_R_IO_H

#include <RcppArmadillo.h>

#include "bayests/priors.h"

#include <string>

// Translation between the R model object and the structs the vendored BayesTS
// core takes. This is the R counterpart of that project's src/io/hdf5/ layer:
// the only place that knows both an `Rcpp::List` and a `bayests::` struct, so
// that the sampler knows neither.
//
// One convention is worth stating once, because getting it wrong is silent
// rather than loud: draws run along the columns inside the core and along the
// rows in R, so everything crossing this boundary goes through draws_to_r() or
// read_draws_if_present().
//
// This is a trimmed copy of the same header in bvartools. What is missing is
// what only a VAR or a VEC needs -- read_spec(), the cointegration space
// priors, the variable selection reader. A dynamic factor model's own reader
// lives next to its binding in DfmNormalGamma.cpp, because its model list names
// its dimensions differently from every other model in that family and the
// translation is not shared with anything.
namespace bayests_r
{

inline bool has(const Rcpp::List &list, const char *name)
{
  return list.containsElementNamed(name) && !Rf_isNull(list[name]);
}

inline int optional_int(const Rcpp::List &list, const char *name, int fallback)
{
  return has(list, name) ? Rcpp::as<int>(list[name]) : fallback;
}

inline void read_mat_if_present(const Rcpp::List &list, const char *name, arma::mat &out)
{
  if (has(list, name)) {
    out = Rcpp::as<arma::mat>(list[name]);
  }
}

inline void read_vec_if_present(const Rcpp::List &list, const char *name, arma::vec &out)
{
  if (has(list, name)) {
    out = Rcpp::as<arma::vec>(list[name]);
  }
}

/// Draws as the core keeps them, from the row-per-draw matrix R keeps.
inline void read_draws_if_present(const Rcpp::List &list, const char *name, arma::mat &out)
{
  if (has(list, name)) {
    out = arma::trans(Rcpp::as<arma::mat>(list[name]));
  }
}

/// Draws as R expects them, from the column-per-draw matrix the core returns.
inline arma::mat draws_to_r(const arma::mat &draws)
{
  return arma::trans(draws);
}

/// The (shape, rate) pair every gamma prior is stored as.
inline bayests::GammaPrior read_gamma_prior(const Rcpp::List &group)
{
  bayests::GammaPrior prior;
  read_vec_if_present(group, "shape", prior.shape);
  read_vec_if_present(group, "rate", prior.rate);
  return prior;
}

} // namespace bayests_r

#endif // DFMTOOLS_BAYESTS_R_IO_H
