// SPDX-License-Identifier: GPL-2.0-or-later

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "dfm_r_translation.h"
#include "bayests/dfm_tvp_gamma.h"

// Dynamic factor model whose loadings and factor transition follow random walks,
// with independent gamma priors on both error precisions. The numerics are the
// vendored BayesTS core; what is here is the translation between the R model
// object and the core's structs. See src/core/VENDORED.md.
//
// DfmNormalGamma.cpp with both coefficient blocks widened from a point to a
// path. Everything that file says about the 'dfmodel' list still applies --
// `model$m` is the core's `k` and `model$n` its `n_factors`, the data live at
// `data$x`, the two error precisions arrive as square matrices and go to the core
// as their diagonals, and the free loadings are ordered differently -- so read
// that one first. Three things are new here, and each is silent if got wrong:
//
//   * `initial$lambda` and `initial$a` are matrices, one column per period,
//     where the constant-coefficient models have a vector. Each comes with the
//     precision of its state innovations (`lambda_sigma_inv`, `a_sigma_inv`) and
//     the state before the sample (`lambda_init`, `a_init`), which is what the
//     prior groups' new `shape`/`rate` pair is over. The names follow bvartools'
//     time-varying VAR and VEC models, where the coefficient prior group gains
//     exactly that pair.
//   * the loading permutation now applies to a path and to two more objects. R
//     keeps the free loadings in `lower.tri()` order and the core row by row, so
//     the path's *rows* are permuted, and so are both the rows and the columns of
//     `lambda_sigma_inv`. Permuting the path's columns instead would compile,
//     run, and reorder the sample.
//   * the two entry points want different slices of the same posterior. A
//     forecast holds the coefficients at their last in-sample period, so
//     `lambda` and `a` are cut to it; the pointwise log likelihood scores every
//     period under its own loadings and takes the whole path. The core refuses
//     the wrong one rather than reading the first period out of it, which is why
//     the cut is made here and not left to the caller.
//
// On the way out, `lambda$sigma` is permuted back into R's ordering. It is one
// number per free loading, so unlike `lambda$coeffs` -- which is the whole M x N
// matrix and needs no permutation -- it would otherwise not line up with the
// prior the caller supplied.

namespace {

using namespace bayests_r;
using dfmtools::lambda_row_major_order;

/// The diagonal of a precision that R stores as a full matrix. Both of this
/// model's error precisions are diagonal by assumption, and the core carries them
/// as vectors for that reason; `add_initial_values.dfmodel()` writes them square.
arma::vec diagonal_of(const Rcpp::List &list, const char *name) {
  arma::mat square;
  read_mat_if_present(list, name, square);
  return square.n_elem == 0 ? arma::vec() : arma::vec(square.diag());
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

bayests::DfmTvpGammaInput read_input(const Rcpp::List &object) {

  bayests::DfmTvpGammaInput input;

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
    if (has(priors, "u")) {
      input.u_sigma_prior = read_gamma_prior(Rcpp::List(priors["u"]));
    }
    if (has(priors, "v")) {
      input.v_sigma_prior = read_gamma_prior(Rcpp::List(priors["v"]));
    }
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

    input.initial.u_sigma_inv = diagonal_of(initial, "uinv");
    input.initial.v_sigma_inv = diagonal_of(initial, "vinv");
  }

  return input;
}

/// Periods in the sample, or zero for an object that carries none. Both slices
/// below are cut from the end of a path, so this is the number they depend on.
int sample_periods(const bayests::DfmTvpGammaInput &input) {
  return input.train.y.n_elem == 0 || input.spec.k <= 0
           ? 0
           : static_cast<int>(input.train.periods(input.spec.k));
}

/// The draws both entry points read, and the one choice between them.
///
/// The factor path is always taken whole: it is part of this posterior rather
/// than derivable from it, so neither entry point works without it, and a
/// forecast starts from the last p factors rather than from one period.
///
/// `terminal` cuts the two coefficient paths to their last in-sample period,
/// which is what a forecast holds them at. The log likelihood passes false and
/// gets the whole path, every period under its own loadings.
bayests::DfmTvpGammaDraws read_draws(const Rcpp::List &object,
                                     const bayests::DfmTvpGammaInput &input,
                                     const bool terminal) {

  bayests::DfmTvpGammaDraws draws;

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
/// coefficient block. `lambda` and `a` are paths here -- m * n * tt and
/// n_a * tt columns per draw against m * n and n_a there -- with the periods
/// stacked within a row.
///
/// `lambda$sigma` is one number per free loading and goes back in R's ordering;
/// `lambda$coeffs` is the whole M x N matrix per period and is already in it.
Rcpp::List write_draws(const bayests::DfmTvpGammaDraws &draws, const arma::uvec &order) {

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

// [[Rcpp::export(.DfmTvpGammaCoefficients)]]
Rcpp::List DfmTvpGammaCoefficients(Rcpp::List object) {

  const bayests::DfmTvpGammaInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  dfmtools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::DfmTvpGammaDraws draws =
    bayests::DfmTvpGammaSampler().draw_coefficients(input, reporter);

  const arma::uvec order = lambda_row_major_order(input.spec.k, input.spec.n_factors);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws, order));
}

// [[Rcpp::export(.DfmTvpGammaForecasts)]]
Rcpp::List DfmTvpGammaForecasts(Rcpp::List object) {

  const bayests::DfmTvpGammaInput input = read_input(object);
  const bayests::DfmTvpGammaDraws draws = read_draws(object, input, true);

  dfmtools::RcppReporter reporter;

  // The loadings and the transition are held at their last in-sample period over
  // the horizon, which is what every time-varying model in this family does and
  // what the posterior supports: the variance of the state innovations is a state
  // of the chain rather than something the draws carry, so there is nothing to
  // extrapolate the random walk with.
  const bayests::ForecastDraws forecast =
    bayests::DfmTvpGammaSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.DfmTvpGammaLogLik)]]
Rcpp::List DfmTvpGammaLogLik(Rcpp::List object) {

  const bayests::DfmTvpGammaInput input = read_input(object);
  const bayests::DfmTvpGammaDraws draws = read_draws(object, input, false);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary. Conditional on the drawn
  // factor path, and scored period by period under that period's loadings.
  const arma::mat loglik = bayests::DfmTvpGammaSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("bem_dfmdata")

spec <- function(n = 1, p = 1, iterations = 400, burnin = 200) {
  object <- create_dfmodel(x = bem_dfmdata, p = p, n = n, tvp = TRUE,
                           iterations = iterations, burnin = burnin)
  object <- add_priors(object,
                       lambda = list(vinv = .01, shape = 3, rate = .01),
                       a = list(vinv = .01, shape = 3, rate = .01),
                       u = list(shape = 5, rate = 4),
                       v = list(shape = 5, rate = 4))
  add_initial_values(object)
}

## Draws run along the rows on this side of the boundary, so a posterior mean is
## a column mean. lambda is m * n * tt wide rather than m * n: a path per draw.
object <- .DfmTvpGammaCoefficients(spec())
vapply(object$posterior, function(x) paste(dim(x$coeffs), collapse = " x "), "")

## The loading of the second series on the factor, period by period.
m <- object$model$m
n <- object$model$n
tt <- nrow(object$data$x)
matrix(colMeans(object$posterior$lambda$coeffs), m * n, tt)[2, ]

*/
