#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/dfm_normal_gamma.h"

#include <algorithm>

// Dynamic factor model with a normal prior on the loadings and on the factor
// transition, and independent gamma priors on both error precisions. The
// numerics are the vendored BayesTS core; what is here is the translation
// between the R model object and the core's structs. See src/core/VENDORED.md.
//
// This is the whole of the posterior simulation this package performs, and it
// replaces the R function dfmpost() that used to do it. Three things about the
// 'dfmodel' list have to be translated, all of them silent if got wrong:
//
//   * the dimensions have different names from the core's. `model$m` is the
//     number of observed series, which the core calls `k`; `model$n` is the
//     number of factors, which the core calls `n_factors` and which every VAR
//     and VEC in the same core uses for its unrestricted deterministic terms.
//   * the data live at `data$x`, where the core expects the response of a model
//     with no regressors at all.
//   * the free loadings are ordered differently. See below.
//
// dfmpost() was correct at one factor and a transition of order at most two.
// Three defects put the limit there, and the BayesTS CHANGELOG lists them; this
// reproduces none of them.

namespace {

using namespace bayests_r;

/// R keeps the free loadings in `lower.tri()` order -- column-major, so every
/// free row of column 0, then of column 1, and so on -- while the core wants
/// them row by row, which is the order its equation-by-equation draw consumes
/// them in and what makes each equation's slice of the prior contiguous.
///
/// This is that permutation: element `i` of the result is the R position of the
/// core's `i`-th free loading, so `core = r.elem(order)` and, for the prior
/// precision, `core = r.submat(order, order)`.
///
/// It matters for the prior, not for the starting value alone. dfmpost() got
/// away with holding one order and slicing the other because `add_priors()`
/// builds that precision as `diag(vinv, n_lambda)`, which reads the same either
/// way; a caller who supplies anything else was served a permuted prior there
/// and gets the right one here.
arma::uvec lambda_row_major_order(const int m, const int n) {

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

/// The diagonal of a precision that R stores as a full matrix. Both of this
/// model's error precisions are diagonal by assumption, and the core carries
/// them as vectors for that reason; `add_initial_values.dfmodel()` writes them
/// out square.
arma::vec diagonal_of(const Rcpp::List &list, const char *name) {
  arma::mat square;
  read_mat_if_present(list, name, square);
  return square.n_elem == 0 ? arma::vec() : arma::vec(square.diag());
}

bayests::DfmNormalGammaInput read_input(const Rcpp::List &object) {

  bayests::DfmNormalGammaInput input;

  const Rcpp::List model = object["model"];
  input.spec.k = Rcpp::as<int>(model["m"]);
  input.spec.n_factors = Rcpp::as<int>(model["n"]);
  input.spec.p = optional_int(model, "p", 0);
  input.spec.iterations = Rcpp::as<int>(model["iterations"]);
  input.spec.burnin = Rcpp::as<int>(model["burnin"]);
  input.spec.h = optional_int(model, "h", 0);

  const int m = input.spec.k;
  const int n = input.spec.n_factors;
  const arma::uvec order = lambda_row_major_order(m, n);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    read_mat_if_present(data, "x", input.train.y);
  }

  // Only draw_coefficients() needs the priors and the initial values, so a
  // missing one is left for validate() to complain about if it turns out to
  // matter.
  if (has(object, "priors")) {
    const Rcpp::List priors = object["priors"];

    if (has(priors, "lambda")) {
      const Rcpp::List prior_lambda = priors["lambda"];
      arma::mat v_inv;
      read_mat_if_present(prior_lambda, "vinv", v_inv);
      if (v_inv.n_rows == order.n_elem) {
        input.lambda_prior.v_inv = v_inv.submat(order, order);
      }
      // add_priors.dfmodel() writes no prior mean for the loadings, so the
      // implied one is zero. Read anyway, for a caller that supplies it.
      arma::vec mu;
      read_vec_if_present(prior_lambda, "mu", mu);
      input.lambda_prior.mu = mu.n_elem == order.n_elem ? arma::vec(mu.elem(order))
                                                        : arma::vec(order.n_elem, arma::fill::zeros);
    }

    if (has(priors, "a")) {
      const Rcpp::List prior_a = priors["a"];
      read_vec_if_present(prior_a, "mu", input.a_prior.mu);
      read_mat_if_present(prior_a, "vinv", input.a_prior.v_inv);
    }

    // "u" and "v" rather than "u_sigma" and "v_sigma": the R names.
    if (has(priors, "u")) {
      input.u_sigma_prior = read_gamma_prior(Rcpp::List(priors["u"]));
    }
    if (has(priors, "v")) {
      input.v_sigma_prior = read_gamma_prior(Rcpp::List(priors["v"]));
    }
  }

  if (has(object, "initial")) {
    const Rcpp::List initial = object["initial"];

    arma::vec lambda;
    read_vec_if_present(initial, "lambda", lambda);
    if (lambda.n_elem == order.n_elem) {
      input.initial.lambda = lambda.elem(order);
    }

    read_vec_if_present(initial, "a", input.initial.a);
    input.initial.u_sigma_inv = diagonal_of(initial, "uinv");
    input.initial.v_sigma_inv = diagonal_of(initial, "vinv");
  }

  return input;
}

/// Both the forecast and the log likelihood read every draw as it stands. The
/// factor path is among them: it is part of this posterior rather than
/// derivable from it, so neither works without it.
bayests::DfmNormalGammaDraws read_draws(const Rcpp::List &object) {

  bayests::DfmNormalGammaDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  if (has(posterior, "lambda")) {
    read_draws_if_present(Rcpp::List(posterior["lambda"]), "coeffs", draws.lambda);
  }
  if (has(posterior, "factors")) {
    read_draws_if_present(Rcpp::List(posterior["factors"]), "coeffs", draws.factors);
  }
  if (has(posterior, "a")) {
    read_draws_if_present(Rcpp::List(posterior["a"]), "coeffs", draws.a);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }
  if (has(posterior, "v_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["v_sigma_inv"]), "coeffs", draws.v_sigma_inv);
  }

  return draws;
}

/// `lambda` comes back as the whole M x N matrix, identifying ones and zeros
/// included, so a caller reshapes rather than having to know which elements were
/// free. The two error blocks come back as precisions, as they do for every other
/// model here; dfmpost() reports variances, so a caller after its output takes
/// the reciprocal.
Rcpp::List write_draws(const bayests::DfmNormalGammaDraws &draws) {

  Rcpp::List posteriors =
    Rcpp::List::create(Rcpp::Named("lambda") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.lambda)),
                       Rcpp::Named("factors") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.factors)),
                       Rcpp::Named("a") = R_NilValue,
                       Rcpp::Named("u_sigma_inv") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv)),
                       Rcpp::Named("v_sigma_inv") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.v_sigma_inv)));

  if (draws.has_a()) {
    posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a));
  }

  return posteriors;
}

} // namespace

// [[Rcpp::export(.DfmNormalGammaCoefficients)]]
Rcpp::List DfmNormalGammaCoefficients(Rcpp::List object) {

  const bayests::DfmNormalGammaInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  dfmtools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming
  // the first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::DfmNormalGammaDraws draws =
    bayests::DfmNormalGammaSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.DfmNormalGammaForecasts)]]
Rcpp::List DfmNormalGammaForecasts(Rcpp::List object) {

  const bayests::DfmNormalGammaInput input = read_input(object);
  const bayests::DfmNormalGammaDraws draws = read_draws(object);

  dfmtools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::DfmNormalGammaSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.DfmNormalGammaLogLik)]]
Rcpp::List DfmNormalGammaLogLik(Rcpp::List object) {

  const bayests::DfmNormalGammaInput input = read_input(object);
  const bayests::DfmNormalGammaDraws draws = read_draws(object);

  // Draws by periods already, which is the orientation R wants and the one
  // WAIC and PSIS-LOO expect; no transpose at this boundary. Conditional on the
  // drawn factor path -- see the header of the sampler.
  const arma::mat loglik = bayests::DfmNormalGammaSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("bem_dfmdata")

spec <- function(n = 1, p = 1, iterations = 400, burnin = 200) {
  object <- create_dfmodel(x = bem_dfmdata, p = p, n = n,
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object,
                       lambda = list(vinv = .01),
                       u = list(shape = 5, rate = 4),
                       a = list(vinv = .01),
                       v = list(shape = 5, rate = 4))
  add_initial_values(object)
}

## Draws run along the rows on this side of the boundary, so a posterior mean is
## a column mean.
object <- .DfmNormalGammaCoefficients(spec())
vapply(object$posterior, function(x) paste(dim(x$coeffs), collapse = " x "), "")

## Two factors and a transition of order two, which is past what dfmpost() could
## do and the reason this exists.
two <- .DfmNormalGammaCoefficients(spec(n = 2, p = 2))
matrix(colMeans(two$posterior$lambda$coeffs), ncol = 2)[1:4, ]

*/
