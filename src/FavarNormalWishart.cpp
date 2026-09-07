#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/favar_normal_wishart.h"

// Factor augmented VAR with a normal prior on the loadings and on the state
// transition, independent gamma priors on the idiosyncratic precisions and a
// Wishart prior on the precision of the state innovations. The numerics are the
// vendored BayesTS core; what is here is the translation between the R model
// object and the core's structs. See src/core/VENDORED.md.
//
// Three things differ from the four dynamic factor bindings beside it, and each
// is a place a silent mistake could live:
//
//   * the state is part data. `model$n` is the number of *unobserved* factors
//     and `model$n_obs` the number of observed ones; the core calls the sum
//     n_state and sizes the transition, the loadings and Q against it. A model
//     object that reports one where the other belongs runs and estimates the
//     wrong thing.
//   * the observed block is data, not a regressor. It arrives at `data$y` and
//     goes into `train.f_obs`, the member the core keeps for it. Putting it in
//     `x` or `z` would read as "a FAVAR is a regression with these on the
//     right", which is the one thing it is not.
//   * the free loadings need no permutation. A dynamic factor model orders them
//     around a unit lower triangle and needs lambda_row_major_order() to get
//     between R's ordering and the core's; here every row after the identifying
//     block is free across its whole width, so the free elements are a plain
//     (k - n) x n_state block and row-major flattening is the whole of it.

namespace {

using namespace bayests_r;

/// The diagonal of a precision that R stores as a full matrix. The idiosyncratic
/// precision is diagonal by assumption -- that assumption is what makes this a
/// factor model -- and the core carries it as a vector for that reason.
arma::vec diagonal_of(const Rcpp::List &list, const char *name) {
  arma::mat square;
  read_mat_if_present(list, name, square);
  return square.n_elem == 0 ? arma::vec() : arma::vec(square.diag());
}

/// The free loadings as the core orders them: row by row, left to right, over
/// the rows after the identifying block.
///
/// R holds them as a (k - n) x n_state matrix, which is the shape a reader wants
/// -- one row per series that carries a loading, one column per state element.
/// The core wants that flattened row-major, which `vectorise(t(L))` is.
arma::vec free_loadings_from_r(const Rcpp::List &list, const char *name) {
  arma::mat block;
  read_mat_if_present(list, name, block);
  if (block.n_elem == 0) {
    return arma::vec();
  }
  return arma::vectorise(arma::trans(block));
}

bayests::FavarNormalWishartInput read_input(const Rcpp::List &object) {

  bayests::FavarNormalWishartInput input;

  const Rcpp::List model = object["model"];
  input.spec.k = Rcpp::as<int>(model["m"]);
  input.spec.n_factors = Rcpp::as<int>(model["n"]);
  input.spec.n_obs_factors = optional_int(model, "n_obs", 0);
  input.spec.p = optional_int(model, "p", 0);
  input.spec.iterations = Rcpp::as<int>(model["iterations"]);
  input.spec.burnin = Rcpp::as<int>(model["burnin"]);
  input.spec.h = optional_int(model, "h", 0);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    read_mat_if_present(data, "x", input.train.y);
    // The observed factors, into the member the core keeps for exactly them.
    // Not `x` and not `z`: those are the regressor layouts, and these are not
    // regressors -- they are the observed half of the state, and they appear on
    // the left of the transition as well as the right.
    read_mat_if_present(data, "y", input.train.f_obs);
  }

  if (has(object, "priors")) {
    const Rcpp::List priors = object["priors"];

    if (has(priors, "lambda")) {
      const Rcpp::List prior_lambda = priors["lambda"];
      read_vec_if_present(prior_lambda, "mu", input.lambda_prior.mu);
      read_mat_if_present(prior_lambda, "vinv", input.lambda_prior.v_inv);
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
      const Rcpp::List prior_v = priors["v"];
      input.v_sigma_prior.df = optional_int(prior_v, "df", 0);
      read_mat_if_present(prior_v, "scale", input.v_sigma_prior.scale);
    }
  }

  if (has(object, "initial")) {
    const Rcpp::List initial = object["initial"];

    input.initial.lambda = free_loadings_from_r(initial, "lambda");
    read_vec_if_present(initial, "a", input.initial.a);
    input.initial.u_sigma_inv = diagonal_of(initial, "uinv");
    // Q is unrestricted, so its precision stays a matrix all the way through --
    // the one error block in this package that is not a diagonal.
    read_mat_if_present(initial, "vinv", input.initial.v_sigma_inv);
  }

  return input;
}

/// The forecast and the log likelihood read every draw as it stands. The factor
/// path is among them: it is part of this posterior rather than derivable from
/// it, so neither works without it.
bayests::FavarNormalWishartDraws read_draws(const Rcpp::List &object) {

  bayests::FavarNormalWishartDraws draws;

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

/// `lambda` comes back as the whole k x n_state matrix, identifying ones and
/// zeros included, so a caller reshapes rather than having to know which
/// elements were free. `factors` is the unobserved block alone -- the observed
/// half of the state is the data the caller already has.
Rcpp::List write_draws(const bayests::FavarNormalWishartDraws &draws) {

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

// [[Rcpp::export(.FavarNormalWishartCoefficients)]]
Rcpp::List FavarNormalWishartCoefficients(Rcpp::List object) {

  const bayests::FavarNormalWishartInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  dfmtools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::FavarNormalWishartDraws draws =
    bayests::FavarNormalWishartSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.FavarNormalWishartForecasts)]]
Rcpp::List FavarNormalWishartForecasts(Rcpp::List object) {

  const bayests::FavarNormalWishartInput input = read_input(object);
  const bayests::FavarNormalWishartDraws draws = read_draws(object);

  dfmtools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::FavarNormalWishartSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.FavarNormalWishartLogLik)]]
Rcpp::List FavarNormalWishartLogLik(Rcpp::List object) {

  const bayests::FavarNormalWishartInput input = read_input(object);
  const bayests::FavarNormalWishartDraws draws = read_draws(object);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary.
  const arma::mat loglik = bayests::FavarNormalWishartSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}
