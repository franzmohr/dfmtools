// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_DFM_SUPPORT_H
#define BAYESTS_CORE_MODELS_DFM_SUPPORT_H

#include "bayests/arma.h"
#include "core/algorithms/chan_jeliazkov_2009.h"
#include "core/models/model_support.h"

#include <algorithm>

namespace bayests::core
{

/// What separates a dynamic factor model from a regression: the regressors are
/// unobserved, so every draw begins by producing a whole tt-period path of them,
/// and the matrix that maps them to the data is identified only up to a rotation
/// and a scale, so part of it is fixed rather than drawn.
///
/// Everything here is one of the two conventions that follow from that pair, and
/// both are easy to get subtly wrong in a way that still runs:
///
///   - Lambda's leading N x N block is unit lower triangular. That leaves
///     min(i, N) free elements in row i, so the rows have *different* widths and
///     the free elements are addressed row by row, not as a rectangle.
///   - Factors before the sample are zero rather than drawn. The transition
///     therefore holds from the first period, in a truncated form for the first
///     p of them, and that truncation is what the path's prior is.
///
/// `factors` is N x tt throughout -- one period per column -- and `a_mat` is the
/// transition blocks side by side, N x Np, so A_j is columns (j-1)N .. jN-1.

/// The observed series as the samplers work with them: M x tt, one period per
/// column.
///
/// Built through stacked_response() rather than by transposing `train.y`,
/// because `y` is allowed to arrive already stacked -- a single row or a single
/// column of vec(y') is how the HDF5 files store it -- and only the stacked
/// vector is the same object in all three cases.
inline arma::mat response_by_period(const TrainData &train, const int k, const int tt)
{
    return arma::reshape(stacked_response(train), k, tt);
}

/// Free elements in row `i` of an M x N loading matrix: everything left of the
/// diagonal while the row is inside the identifying block, the whole row after
/// it. Row zero has none -- its only entry is the leading one.
inline int lambda_row_width(const int i, const int n_factors)
{
    return std::min(i, n_factors);
}

/// An M x N loading matrix carrying nothing but its identification: ones on the
/// diagonal of the leading block, zeros everywhere else.
inline arma::mat identified_loadings(const int k, const int n_factors)
{
    arma::mat lambda(k, n_factors, arma::fill::zeros);
    lambda.diag().ones();
    return lambda;
}

/// Writes the free loadings into Lambda, leaving the identifying block alone.
///
/// The vector runs row by row and, within a row, left to right -- the order the
/// equation-by-equation draw consumes it in, which is what makes each row's
/// slice of the prior contiguous. Not the order `lower.tri()` produces in R,
/// which is column-major; the two coincide only at one factor, and bvartools'
/// dfmpost() gets away with the mismatch because its prior precision is a scalar
/// diagonal and so reads the same either way.
inline void fill_lambda(arma::mat &lambda, const arma::vec &free)
{
    const int k = static_cast<int>(lambda.n_rows);
    const int n = static_cast<int>(lambda.n_cols);
    int pos = 0;
    for (int i = 1; i < k; i++)
    {
        const int width = lambda_row_width(i, n);
        lambda.submat(i, 0, i, width - 1) = arma::trans(free.subvec(pos, pos + width - 1));
        pos += width;
    }
}

/// The lagged factors the transition regresses on: (N p) x tt, block j holding
/// f_{t-j}, with the columns before the sample left at zero.
///
/// Block j occupies rows (j-1)N .. jN-1. bvartools' dfmpost() writes rows
/// (j-1) .. (j-1)+N-1 instead, which is the same thing only at one factor and
/// otherwise lays the blocks on top of one another.
inline void fill_lagged_factors(arma::mat &x_a, const arma::mat &factors, const int n, const int p)
{
    const int tt = static_cast<int>(factors.n_cols);
    x_a.zeros();
    for (int j = 1; j <= p; j++)
    {
        x_a.submat((j - 1) * n, j, j * n - 1, tt - 1) = factors.cols(0, tt - 1 - j);
    }
}

/// The transition residuals f_t - sum_j A_j f_{t-j}, N x tt, on the same
/// zero-before-the-sample convention.
inline arma::mat transition_residuals(const arma::mat &factors, const arma::mat &a_mat,
                                      const int n, const int p)
{
    const int tt = static_cast<int>(factors.n_cols);
    arma::mat v = factors;
    for (int j = 1; j <= p; j++)
    {
        v.cols(j, tt - 1) -= a_mat.cols((j - 1) * n, j * n - 1) * factors.cols(0, tt - 1 - j);
    }
    return v;
}

/// One draw of a diagonal precision under independent gamma priors: for each row
/// of `residuals`, Gamma(shape + T/2, rate + sum_t e_t^2 / 2). `post_shape` is
/// the first of those, which does not change over the chain and is formed once.
///
/// Row-wise rather than off the diagonal of `residuals * residuals'`, which is
/// the spelling the constant-variance VARs use. There k is small; here it is the
/// number of observed series and a factor model is worth having precisely when
/// that is large, so forming the k x k product to read k numbers off it would be
/// the most expensive thing in the block.
inline void draw_diagonal_precision(arma::vec &precision, const arma::mat &residuals,
                                    const arma::vec &post_shape, const arma::vec &prior_rate)
{
    const arma::vec sse = arma::sum(arma::square(residuals), 1);
    for (arma::uword i = 0; i < precision.n_elem; i++)
    {
        precision(i) = arma::randg<double>(
            arma::distr_param(post_shape(i), 1.0 / (prior_rate(i) + sse(i) * 0.5)));
    }
}

/// Adds what the data contribute to the factor transition's posterior when the
/// factor innovations' covariance moves with time. `precision` and `rhs` come in
/// seeded with the prior's contribution and go out with both.
///
/// The constant-covariance case collapses: sum_t kron(x_t x_t', V^-1) is
/// kron(X X', V^-1), which is why DfmNormalGamma never forms a per-period
/// anything. With a V_t per period it does not, and the sum has to be taken.
///
/// Taking it as a sum of tt Kronecker products would be the direct spelling and
/// is what the test checks this against. It is not what happens here, because V_t
/// is diagonal: the N equations of the transition are then conditionally
/// independent, so equation i's own weighted cross-product X W_i X' -- where W_i
/// is the diagonal of period-by-period precisions of factor i -- is the whole of
/// what the data give it, and the off-diagonal blocks that a Kronecker product
/// would spend time filling with zeros are zero. That turns tt products of size
/// (Np)^2 N^2 into N products of size (Np)^2 tt, and leaves nothing allocated per
/// period.
///
/// The scatter is the one piece of index arithmetic worth reading twice.
/// `a` is vec([A_1 .. A_p]) of an N x Np matrix, column-major, so element (i, c)
/// of that matrix sits at c N + i -- row i of A occupies positions i, N + i,
/// 2N + i, and so on. Entry (c, d) of equation i's cross-product therefore
/// belongs at (c N + i, d N + i).
///
/// @param precision (N^2 p) square, added to in place.
/// @param rhs (N^2 p), added to in place.
/// @param x_a (N p) x tt lagged factors, as fill_lagged_factors() writes them.
/// @param factors N x tt factor path.
/// @param v_precision tt x N, one row per period: the diagonal of V_t^-1.
inline void accumulate_transition_moments(arma::mat &precision, arma::vec &rhs,
                                          const arma::mat &x_a, const arma::mat &factors,
                                          const arma::mat &v_precision, const int n, const int p)
{
    const int np = n * p;

    for (int i = 0; i < n; i++)
    {
        arma::mat x_weighted = x_a;
        x_weighted.each_row() %= arma::trans(v_precision.col(i));

        const arma::mat cross = x_weighted * arma::trans(x_a);
        const arma::vec moment = x_weighted * arma::trans(factors.row(i));

        for (int c = 0; c < np; c++)
        {
            const int row = c * n + i;
            rhs(row) += moment(c);
            for (int d = 0; d < np; d++)
            {
                precision(row, d * n + i) += cross(c, d);
            }
        }
    }
}

/// Covariance of the first p factors, which is what the path's prior is.
///
/// The factors before the sample are zero rather than drawn, so the first p
/// periods are not a free prior at all -- they are the transition itself, run
/// from nothing:
///
///     f_1 = v_1,   f_2 = A_1 f_1 + v_2,   ...,   f_p = sum_{j<p} A_j f_{p-j} + v_p.
///
/// Stacked, that is H f = v with H unit block lower triangular carrying -A_j on
/// its j-th subdiagonal, so Cov(f_{1..p}) = H^-1 Cov(v_{1..p}) H^-T. Handing that
/// to chan_jeliazkov_2009 as `P_init` is what makes its prior-plus-transitions
/// decomposition reproduce this model exactly rather than approximately: the
/// blocks it then builds are the same ones the (tt N) square precision of the
/// reference implementation has.
///
/// H is unit triangular, hence exactly invertible, and every V is positive
/// definite, so the result is too -- which chan_jeliazkov_2009 requires of it.
///
/// `v_sigma` is either one N x N covariance that holds in every period, giving
/// the block diagonal I_p kron V, or a (p N) x N stack of the first p periods'
/// own covariances. The second is what a model whose factor innovations carry
/// stochastic volatility has, and taking the first period's for all p of them
/// there would misstate the prior by however much the volatility moved over the
/// first p periods -- silently, since either shape runs. Dispatching on the
/// height is the same arrangement chan_jeliazkov_2009 uses for the same reason.
inline arma::mat initial_state_covariance(const arma::mat &a_mat, const arma::mat &v_sigma,
                                          const int n, const int p)
{
    const int side = p * n;
    arma::mat h = arma::eye<arma::mat>(side, side);
    for (int i = 1; i < p; i++)
    {
        for (int j = 0; j < i; j++)
        {
            h.submat(i * n, j * n, (i + 1) * n - 1, (j + 1) * n - 1) =
                -a_mat.cols((i - j - 1) * n, (i - j) * n - 1);
        }
    }

    const int v_stride = (static_cast<int>(v_sigma.n_rows) == n) ? 0 : n;
    arma::mat v_blocks(side, side, arma::fill::zeros);
    for (int i = 0; i < p; i++)
    {
        v_blocks.submat(i * n, i * n, (i + 1) * n - 1, (i + 1) * n - 1) =
            v_sigma.rows(i * v_stride, i * v_stride + n - 1);
    }

    const arma::mat h_inv = arma::solve(arma::trimatl(h), arma::eye<arma::mat>(side, side));
    return arma::symmatu(h_inv * v_blocks * arma::trans(h_inv));
}

/// One draw of the whole factor path, N x tt.
///
/// The path is a single Gaussian vector whose precision is block banded of
/// bandwidth p -- a period's measurement touches one factor, and the transition
/// couples p of them -- so chan_jeliazkov_2009 draws it in one sweep over the
/// periods. The alternative, which bvartools' dfmpost() takes, is to form the
/// (tt N) square precision and factorise it: the same distribution at
/// O(tt^3 N^3) instead of O(tt N^3).
///
/// The loading matrix is the measurement design and is the same in every period,
/// which is the case chan_jeliazkov_2009 assembles once rather than tt times --
/// worth having here, since a factor model is worth having when M is large.
///
/// `a_mat` is N x Np for a transition of order p >= 1. A model whose factors have
/// no dynamics passes an N x N block of zeros: f_t = 0 f_{t-1} + v_t *is* the
/// serially independent factor it then has, and its prior on the first state is
/// V, so the two cases are one code path rather than a block-diagonal special
/// case that would draw the same thing.
///
/// The last of the tt + 1 columns chan_jeliazkov_2009 returns is dropped. It is
/// the transition applied once past the end of the sample and no observation
/// informs it, so integrating it out leaves the rest exactly as it was.
inline arma::mat draw_factor_path(const arma::mat &x_t, const arma::mat &lambda,
                                  const arma::mat &u_sigma, const arma::mat &v_sigma,
                                  const arma::mat &a_mat, const int n, const int p_state)
{
    const arma::uword tt = x_t.n_cols;
    const arma::vec a_init(static_cast<arma::uword>(p_state * n), arma::fill::zeros);

    return chan_jeliazkov_2009(x_t, lambda, u_sigma, v_sigma, a_mat, a_init,
                               initial_state_covariance(a_mat, v_sigma, n, p_state))
        .cols(0, tt - 1);
}

/// Writes the per-period diagonal covariances of an error term into the stack
/// chan_jeliazkov_2009 reads: one K x K block per period, block t holding
/// diag(variance.row(t)).
///
/// `stack` must be (K tt) x K and must already be zero. Nothing here writes an
/// off-diagonal element, so a caller that zeroes the buffer once before the
/// chain starts pays tt K stores per draw instead of tt K^2 -- at 100 series over
/// 300 periods, thirty thousand against three million, for a matrix whose
/// off-diagonal was zero throughout.
///
/// `variance` is tt x K, one period per row: the layout the stochastic
/// volatility block works in, so `arma::exp(h)` goes straight in.
inline void fill_stacked_diagonal(arma::mat &stack, const arma::mat &variance)
{
    const arma::uword tt = variance.n_rows;
    const arma::uword k = variance.n_cols;
    for (arma::uword t = 0; t < tt; t++)
    {
        const arma::uword base = t * k;
        for (arma::uword i = 0; i < k; i++)
        {
            stack(base + i, i) = variance(t, i);
        }
    }
}

/// One draw of the whole factor path when both error terms carry stochastic
/// volatility, N x tt.
///
/// `u_stack` is (K tt) x K and `v_stack` is (N tt) x N, both in this model's
/// indexing: block t is period t's covariance. `fill_stacked_diagonal` builds
/// them. Everything draw_factor_path() says still applies -- the path is one
/// Gaussian vector of block banded precision, drawn in a sweep over the periods
/// -- and the loading matrix is still the same in every period, so the
/// measurement design is still assembled once.
///
/// What is different is the indexing, and it is the one thing here that is easy
/// to get wrong in a way that still runs. chan_jeliazkov_2009 indexes the
/// transition that *produces* state column t by t - 1; unit_chan_jeliazkov.cpp
/// pins that convention and kalman_durbin_koopman_2002 shares it. This model's
/// column t is period t and its innovation is v_t ~ N(0, V_t), so block s of the
/// stack handed over must hold V_{s+1} -- `v_stack` moved up by one period.
/// Handing `v_stack` over unshifted would estimate a model whose volatility lags
/// its own innovations by a period, which is a different model and not a broken
/// one: nothing would fail.
///
/// The shift leaves the last block wanting V_tt, which does not exist. It
/// belongs to state column tt, the one past the end of the sample that is
/// dropped on the way out, and any positive definite matrix does there: f_tt
/// enters the joint through that single transition and nothing else, so
/// integrating it out is a Gaussian integral in f_tt alone and leaves the
/// columns that are kept exactly as they were, whatever its covariance. The
/// previous period's is reused because it is already the right shape and already
/// positive definite.
///
/// The prior over the first p columns takes `v_stack` unshifted: those columns
/// are the truncated transitions run from nothing, so the covariances they need
/// are V_0 ... V_{p-1}, which is where the stack starts.
inline arma::mat draw_factor_path_sv(const arma::mat &x_t, const arma::mat &lambda,
                                     const arma::mat &u_stack, const arma::mat &v_stack,
                                     const arma::mat &a_mat, const int n, const int p_state)
{
    const arma::uword tt = x_t.n_cols;
    const arma::uword nn = static_cast<arma::uword>(n);
    const arma::vec a_init(static_cast<arma::uword>(p_state) * nn, arma::fill::zeros);

    arma::mat v_transitions(v_stack.n_rows, v_stack.n_cols);
    v_transitions.head_rows((tt - 1) * nn) = v_stack.tail_rows((tt - 1) * nn);
    v_transitions.tail_rows(nn) = v_stack.tail_rows(nn);

    const arma::mat p_init = initial_state_covariance(
        a_mat, v_stack.head_rows(static_cast<arma::uword>(p_state) * nn), n, p_state);

    return chan_jeliazkov_2009(x_t, lambda, u_stack, v_transitions, a_mat, a_init, p_init)
        .cols(0, tt - 1);
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_DFM_SUPPORT_H
