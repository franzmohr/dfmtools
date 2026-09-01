#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "dfm_r_translation.h"
#include "bayests/dfm_normal_stochvol.h"

// Dynamic factor model with a normal prior on the loadings and on the factor
// transition, and stochastic volatility in both error terms. The numerics are
// the vendored BayesTS core; what is here is the translation between the R model
// object and the core's structs. See src/core/VENDORED.md.
//
// DfmNormalGamma.cpp with the two gamma priors replaced by two stochastic
// volatility blocks. Everything that file says about the 'dfmodel' list still
// applies -- `model$m` is the core's `k` and `model$n` its `n_factors`, the data
// live at `data$x`, and the free loadings are ordered differently -- so read that
// one first. What is new here is the pair of error blocks, and one thing about
// them is worth stating because it is silent if got wrong:
//
//   * they have different widths. `priors$u` is `m` wide, one log-volatility per
//     observed series; `priors$v` is `n` wide, one per factor. Swapping them
//     gives two well-formed lists of the wrong length, which the core's
//     validate() rejects by name -- but only because it is told which is which
//     here.
//
// The field names follow bvartools' stochastic volatility priors exactly --
// `offset`, `shape`, `rate`, `mu`, `v_inv`, `sigma` -- so that a reader who knows
// the VAR side of the family does not have to learn a second spelling. `sigma`
// among them is the *starting value* of the log-volatility innovation variance
// rather than a prior; it sits in the prior group because that is where
// bvartools' model objects put it.

namespace {

using namespace bayests_r;
using dfmtools::lambda_row_major_order;

/// One of the two stochastic volatility prior groups. They differ in their name
/// and their width and in nothing else, so reading them through one function is
/// what keeps the two from drifting apart.
void read_stochvol_group(const Rcpp::List &priors, const char *name,
                         bayests::StochvolPrior &prior, arma::vec &initial_h_sigma) {

  if (!has(priors, name)) {
    return;
  }

  const Rcpp::List group = priors[name];
  read_vec_if_present(group, "offset", prior.offset);
  read_vec_if_present(group, "shape", prior.state.sigma.shape);
  read_vec_if_present(group, "rate", prior.state.sigma.rate);
  read_vec_if_present(group, "mu", prior.state.initial_state.mu);
  read_mat_if_present(group, "v_inv", prior.state.initial_state.v_inv);

  // A state the sampler redraws every iteration, even though R keeps it next to
  // the prior it is drawn under.
  read_vec_if_present(group, "sigma", initial_h_sigma);
}

bayests::DfmNormalStochvolInput read_input(const Rcpp::List &object) {

  bayests::DfmNormalStochvolInput input;

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

    // "u" and "v" rather than "u_sigma" and "v_sigma": the R names. The first is
    // m wide and the second n.
    read_stochvol_group(priors, "u", input.u_sigma_prior, input.initial.u_h_sigma);
    read_stochvol_group(priors, "v", input.v_sigma_prior, input.initial.v_h_sigma);
  }

  if (has(object, "initial")) {
    const Rcpp::List initial = object["initial"];

    arma::vec lambda;
    read_vec_if_present(initial, "lambda", lambda);
    if (lambda.n_elem == order.n_elem) {
      input.initial.lambda = lambda.elem(order);
    }

    read_vec_if_present(initial, "a", input.initial.a);

    // Named apart rather than sharing bvartools' `h` and `h_init`, because this
    // model has two log-volatility paths and neither is the obvious default.
    read_mat_if_present(initial, "u_h", input.initial.u_h);
    read_vec_if_present(initial, "u_h_init", input.initial.u_h_init);
    read_mat_if_present(initial, "v_h", input.initial.v_h);
    read_vec_if_present(initial, "v_h_init", input.initial.v_h_init);
  }

  return input;
}

/// Both the forecast and the log likelihood read every draw as it stands. The
/// factor path is among them: it is part of this posterior rather than derivable
/// from it, so neither works without it.
bayests::DfmNormalStochvolDraws read_draws(const Rcpp::List &object) {

  bayests::DfmNormalStochvolDraws draws;

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

/// The same five elements DfmNormalGamma returns, and the two error blocks under
/// the same names -- but a whole path per draw rather than one number per series,
/// so `u_sigma_inv` is m * tt columns wide here and m there. Precisions, as
/// everywhere in this family.
Rcpp::List write_draws(const bayests::DfmNormalStochvolDraws &draws) {

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

// [[Rcpp::export(.DfmNormalStochvolCoefficients)]]
Rcpp::List DfmNormalStochvolCoefficients(Rcpp::List object) {

  const bayests::DfmNormalStochvolInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  dfmtools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::DfmNormalStochvolDraws draws =
    bayests::DfmNormalStochvolSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.DfmNormalStochvolForecasts)]]
Rcpp::List DfmNormalStochvolForecasts(Rcpp::List object) {

  const bayests::DfmNormalStochvolInput input = read_input(object);
  const bayests::DfmNormalStochvolDraws draws = read_draws(object);

  dfmtools::RcppReporter reporter;

  // Both volatilities are held at their last in-sample value over the horizon,
  // which is what every stochastic volatility model in this family does and what
  // the posterior supports: the variance of the log-volatility innovations is a
  // state of the chain rather than something the draws carry, so there is nothing
  // to extrapolate the random walk with.
  const bayests::ForecastDraws forecast =
    bayests::DfmNormalStochvolSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.DfmNormalStochvolLogLik)]]
Rcpp::List DfmNormalStochvolLogLik(Rcpp::List object) {

  const bayests::DfmNormalStochvolInput input = read_input(object);
  const bayests::DfmNormalStochvolDraws draws = read_draws(object);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary. Conditional on the drawn
  // factor path *and* on the drawn volatility -- see the header of the sampler.
  const arma::mat loglik = bayests::DfmNormalStochvolSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("bem_dfmdata")

spec <- function(n = 1, p = 1, iterations = 400, burnin = 200) {
  object <- create_dfmodel(x = bem_dfmdata, p = p, n = n, error = "sv",
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object,
                       lambda = list(vinv = .01),
                       a = list(vinv = .01),
                       u = list(mu = 0, v_i = 0.1, shape = 3, rate = 0.2,
                                state_variance = 0.05, offset = 1e-4),
                       v = list(mu = 0, v_i = 0.1, shape = 3, rate = 0.2,
                                state_variance = 0.05, offset = 1e-4))
  add_initial_values(object)
}

## Draws run along the rows on this side of the boundary, so a posterior mean is
## a column mean. u_sigma_inv is m * tt wide rather than m: a path per draw.
object <- .DfmNormalStochvolCoefficients(spec())
vapply(object$posterior, function(x) paste(dim(x$coeffs), collapse = " x "), "")

## The idiosyncratic variance of the first series, period by period.
m <- object$model$m
tt <- nrow(object$data$x)
matrix(1 / colMeans(object$posterior$u_sigma_inv$coeffs), m, tt)[1, ]

*/
