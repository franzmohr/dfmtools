// SPDX-License-Identifier: GPL-2.0-or-later

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "dfm_r_translation.h"
#include "bayests/dfm_tvp_stochvol.h"

// Dynamic factor model whose loadings, factor transition and both error
// covariances all move with time. The numerics are the vendored BayesTS core;
// what is here is the translation between the R model object and the core's
// structs. See src/core/VENDORED.md.
//
// DfmTvpGamma.cpp's coefficient blocks over DfmNormalStochvol.cpp's error
// blocks, and nothing that is not in one of those two. Read DfmTvpGamma.cpp for
// the coefficient paths -- the matrices where the constant models have vectors,
// the permutation that then applies to a path's *rows* and to both dimensions of
// `lambda_sigma_inv`, and the different slices the two entry points want -- and
// DfmNormalStochvol.cpp for the two volatility groups and the widths that keep
// them apart.
//
// What is worth stating once here, because it is only visible when the two are
// together: `shape` and `rate` now appear in all four prior groups and mean the
// same thing in each -- the inverse gamma on a random walk's innovation variance
// -- at four different widths. `lambda` is n_lambda wide, `a` is n_factor_a wide,
// `u` is m wide and `v` is n wide. Fetching one group's pair for another gives a
// well-formed list of the wrong length, which the core's validate() rejects by
// name, but only because it is told which is which here.
//
// The log likelihood takes one more thing whole than its gamma sibling does. Both
// want the loading path, because every period is scored under its own Lambda_t;
// this one also wants the idiosyncratic precision path, because every period is
// scored under its own U_t. The forecast takes neither: it holds all four at
// their last in-sample period, and the two precisions need no cut because the
// sampler takes the last block of whatever it is handed.

namespace {

using namespace bayests_r;
using dfmtools::lambda_row_major_order;

/// One of the two stochastic volatility prior groups. They differ in their name
/// and their width and in nothing else, so reading them through one function is
/// what keeps the two from drifting apart. `sigma` among them is the *starting
/// value* of the log-volatility innovation variance rather than a prior; it sits
/// in the prior group because that is where bvartools' model objects put it.
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

/// The state equation of one coefficient block: an inverse gamma on the variance
/// of the state innovations, and a normal on the state of the period before the
/// sample. One function for both blocks, so that the loadings and the transition
/// cannot end up reading a field the other has not got.
///
/// `vinv` rather than `v_inv`, which is what this package calls the prior
/// precision of these two groups and what the constant-coefficient bindings
/// already read.
void read_random_walk_prior(const Rcpp::List &group, bayests::RandomWalkPrior &prior) {

  read_vec_if_present(group, "shape", prior.sigma.shape);
  read_vec_if_present(group, "rate", prior.sigma.rate);
  read_vec_if_present(group, "mu", prior.initial_state.mu);
  read_mat_if_present(group, "vinv", prior.initial_state.v_inv);
}

bayests::DfmTvpStochvolInput read_input(const Rcpp::List &object) {

  bayests::DfmTvpStochvolInput input;

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
      read_random_walk_prior(prior_lambda, input.lambda_prior);

      // The permutation, applied to whatever came back. The prior mean is
      // implied zero where add_priors.dfmodel() wrote none.
      if (input.lambda_prior.initial_state.v_inv.n_rows == order.n_elem) {
        input.lambda_prior.initial_state.v_inv =
          input.lambda_prior.initial_state.v_inv.submat(order, order);
      }
      input.lambda_prior.initial_state.mu =
        input.lambda_prior.initial_state.mu.n_elem == order.n_elem
          ? arma::vec(input.lambda_prior.initial_state.mu.elem(order))
          : arma::vec(order.n_elem, arma::fill::zeros);
      if (input.lambda_prior.sigma.shape.n_elem == order.n_elem) {
        input.lambda_prior.sigma.shape = arma::vec(input.lambda_prior.sigma.shape.elem(order));
        input.lambda_prior.sigma.rate = arma::vec(input.lambda_prior.sigma.rate.elem(order));
      }
    }

    if (has(priors, "a")) {
      read_random_walk_prior(Rcpp::List(priors["a"]), input.a_prior);
    }

    // "u" and "v" rather than "u_sigma" and "v_sigma": the R names. The first is
    // m wide and the second n.
    read_stochvol_group(priors, "u", input.u_sigma_prior, input.initial.u_h_sigma);
    read_stochvol_group(priors, "v", input.v_sigma_prior, input.initial.v_h_sigma);
  }

  if (has(object, "initial")) {
    const Rcpp::List initial = object["initial"];

    // The loading path, one period per column. Rows are the free loadings, so
    // the permutation acts on rows alone -- the columns are the sample.
    arma::mat lambda;
    read_mat_if_present(initial, "lambda", lambda);
    if (lambda.n_rows == order.n_elem) {
      input.initial.lambda = lambda.rows(order);
    }

    arma::mat lambda_sigma_inv;
    read_mat_if_present(initial, "lambda_sigma_inv", lambda_sigma_inv);
    if (lambda_sigma_inv.n_rows == order.n_elem) {
      input.initial.lambda_sigma_inv = lambda_sigma_inv.submat(order, order);
    }

    arma::vec lambda_init;
    read_vec_if_present(initial, "lambda_init", lambda_init);
    if (lambda_init.n_elem == order.n_elem) {
      input.initial.lambda_init = lambda_init.elem(order);
    }

    // The transition path needs no permutation: vec([A_1 .. A_p]) is the same
    // object on both sides.
    read_mat_if_present(initial, "a", input.initial.a);
    read_mat_if_present(initial, "a_sigma_inv", input.initial.a_sigma_inv);
    read_vec_if_present(initial, "a_init", input.initial.a_init);

    // Named apart rather than sharing bvartools' `h` and `h_init`, because this
    // model has two log-volatility paths and neither is the obvious default.
    read_mat_if_present(initial, "u_h", input.initial.u_h);
    read_vec_if_present(initial, "u_h_init", input.initial.u_h_init);
    read_mat_if_present(initial, "v_h", input.initial.v_h);
    read_vec_if_present(initial, "v_h_init", input.initial.v_h_init);
  }

  return input;
}

/// Periods in the sample, or zero for an object that carries none. Both slices
/// below are cut from the end of a path, so this is the number they depend on.
int sample_periods(const bayests::DfmTvpStochvolInput &input) {
  return input.train.y.n_elem == 0 || input.spec.k <= 0
           ? 0
           : static_cast<int>(input.train.periods(input.spec.k));
}

/// The draws both entry points read, and the one choice between them.
///
/// The factor path and the two precision paths are always taken whole. The
/// factors are part of this posterior rather than derivable from it, so neither
/// entry point works without them and a forecast starts from the last p of them;
/// the precisions are taken whole because the sampler reads the last block of
/// whatever it is handed, so a forecast needs no cut and the log likelihood needs
/// the rest.
///
/// `terminal` cuts the two *coefficient* paths to their last in-sample period,
/// which is what a forecast holds them at. The log likelihood passes false and
/// gets the whole path, every period under its own loadings.
bayests::DfmTvpStochvolDraws read_draws(const Rcpp::List &object,
                                     const bayests::DfmTvpStochvolInput &input,
                                     const bool terminal) {

  bayests::DfmTvpStochvolDraws draws;

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

  if (!terminal) {
    return draws;
  }

  const int tt = sample_periods(input);
  const arma::uword lambda_width =
    static_cast<arma::uword>(input.spec.k) * input.spec.n_factors;
  const arma::uword a_width = static_cast<arma::uword>(input.spec.n_factor_a());

  // Cut only what is actually a path. A caller who has already cut it, or a
  // posterior that came from somewhere else, is left alone and the core's own
  // check on the width is what speaks.
  if (tt > 0 && draws.lambda.n_rows == lambda_width * tt) {
    draws.lambda = draws.lambda.tail_rows(lambda_width);
  }
  if (tt > 0 && a_width > 0 && draws.a.n_rows == a_width * tt) {
    draws.a = draws.a.tail_rows(a_width);
  }

  return draws;
}

/// The five elements DfmNormalGamma returns, plus the state variance of each
/// coefficient block. Every one of them but the two state variances is a path:
/// `lambda` is m * n * tt columns per draw, `a` is n_a * tt, `u_sigma_inv` is
/// m * tt and `v_sigma_inv` is n * tt, with the periods stacked within a row.
///
/// `lambda$sigma` is one number per free loading and goes back in R's ordering;
/// `lambda$coeffs` is the whole M x N matrix per period and is already in it.
Rcpp::List write_draws(const bayests::DfmTvpStochvolDraws &draws, const arma::uvec &order) {

  arma::mat lambda_sigma(draws.lambda_sigma.n_rows, draws.lambda_sigma.n_cols);
  if (draws.lambda_sigma.n_rows == order.n_elem) {
    lambda_sigma.rows(order) = draws.lambda_sigma;
  } else {
    lambda_sigma = draws.lambda_sigma;
  }

  Rcpp::List posteriors =
    Rcpp::List::create(Rcpp::Named("lambda") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.lambda),
                         Rcpp::Named("sigma") = draws_to_r(lambda_sigma)),
                       Rcpp::Named("factors") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.factors)),
                       Rcpp::Named("a") = R_NilValue,
                       Rcpp::Named("u_sigma_inv") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv)),
                       Rcpp::Named("v_sigma_inv") = Rcpp::List::create(
                         Rcpp::Named("coeffs") = draws_to_r(draws.v_sigma_inv)));

  if (draws.has_a()) {
    posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a),
                                         Rcpp::Named("sigma") = draws_to_r(draws.a_sigma));
  }

  return posteriors;
}

} // namespace

// [[Rcpp::export(.DfmTvpStochvolCoefficients)]]
Rcpp::List DfmTvpStochvolCoefficients(Rcpp::List object) {

  const bayests::DfmTvpStochvolInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  dfmtools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::DfmTvpStochvolDraws draws =
    bayests::DfmTvpStochvolSampler().draw_coefficients(input, reporter);

  const arma::uvec order = lambda_row_major_order(input.spec.k, input.spec.n_factors);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws, order));
}

// [[Rcpp::export(.DfmTvpStochvolForecasts)]]
Rcpp::List DfmTvpStochvolForecasts(Rcpp::List object) {

  const bayests::DfmTvpStochvolInput input = read_input(object);
  const bayests::DfmTvpStochvolDraws draws = read_draws(object, input, true);

  dfmtools::RcppReporter reporter;

  // The loadings and the transition are held at their last in-sample period over
  // the horizon, which is what every time-varying model in this family does and
  // what the posterior supports: the variance of the state innovations is a state
  // of the chain rather than something the draws carry, so there is nothing to
  // extrapolate the random walk with.
  const bayests::ForecastDraws forecast =
    bayests::DfmTvpStochvolSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.DfmTvpStochvolLogLik)]]
Rcpp::List DfmTvpStochvolLogLik(Rcpp::List object) {

  const bayests::DfmTvpStochvolInput input = read_input(object);
  const bayests::DfmTvpStochvolDraws draws = read_draws(object, input, false);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary. Conditional on the drawn
  // factor path, and scored period by period under that period's loadings.
  const arma::mat loglik = bayests::DfmTvpStochvolSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("bem_dfmdata")

spec <- function(n = 1, p = 1, iterations = 400, burnin = 200) {
  object <- create_dfmodel(x = bem_dfmdata, p = p, n = n, error = "sv", tvp = TRUE,
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object,
                       lambda = list(vinv = .01, shape = 3, rate = .01),
                       a = list(vinv = .01, shape = 3, rate = .01),
                       u = list(mu = 0, v_i = 0.1, shape = 3, rate = 0.2,
                                state_variance = 0.05, offset = 1e-4),
                       v = list(mu = 0, v_i = 0.1, shape = 3, rate = 0.2,
                                state_variance = 0.05, offset = 1e-4))
  add_initial_values(object)
}

## Draws run along the rows on this side of the boundary, so a posterior mean is
## a column mean. Everything but the two state variances is a path per draw.
object <- .DfmTvpStochvolCoefficients(spec())
vapply(object$posterior, function(x) paste(dim(x$coeffs), collapse = " x "), "")

m <- object$model$m
n <- object$model$n
tt <- nrow(object$data$x)

## The loading of the second series on the factor, period by period, beside its
## idiosyncratic variance over the same periods -- the two this model separates.
matrix(colMeans(object$posterior$lambda$coeffs), m * n, tt)[2, ]
matrix(1 / colMeans(object$posterior$u_sigma_inv$coeffs), m, tt)[2, ]

*/
