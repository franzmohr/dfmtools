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
///
/// `a_mat` dispatches on its height the same way, and for the same reason: one
/// N x Np transition that holds in every period, or a (p N) x Np stack of the
/// first p periods' own. A model whose transition drifts needs the second --
/// f_i is produced by A_i, not by A_1 -- and both shapes run, so the stack is
/// passed unshifted here even though chan_jeliazkov_2009 wants it shifted. These
/// p blocks are the transitions the prior *is*, not transitions producing a
/// later column.
inline arma::mat initial_state_covariance(const arma::mat &a_mat, const arma::mat &v_sigma,
                                          const int n, const int p)
{
    const int side = p * n;
    const int a_stride = (static_cast<int>(a_mat.n_rows) == n) ? 0 : n;
    arma::mat h = arma::eye<arma::mat>(side, side);
    for (int i = 1; i < p; i++)
    {
        for (int j = 0; j < i; j++)
        {
            h.submat(i * n, j * n, (i + 1) * n - 1, (j + 1) * n - 1) =
                -a_mat.submat(i * a_stride, (i - j - 1) * n, i * a_stride + n - 1,
                              (i - j) * n - 1);
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

/// A stack moved up by one period, its last block repeated.
///
/// The one index conversion the factor path needs, and the one that is easiest to
/// get wrong in a way that still runs. chan_jeliazkov_2009 indexes the transition
/// that *produces* state column t by t - 1; every model here indexes block t at
/// period t. So a per-period argument that describes a transition -- the
/// coefficients, or the covariance of its innovation -- has to move up one block
/// on the way in. Handing one over unshifted estimates a model whose transitions
/// lag their own period by one, which is a different model and not a broken one:
/// nothing would fail.
///
/// The last block is left wanting the period after the sample, which does not
/// exist. It belongs to state column tt, the one past the end that is dropped on
/// the way out, and anything admissible does there: f_tt enters the joint through
/// that single Gaussian factor and nothing else, so integrating it out is an
/// integral in f_tt alone and leaves the kept columns exactly as they were. The
/// previous period's block is reused because it is already the right shape and
/// already admissible.
inline arma::mat shifted_by_one_period(const arma::mat &stack, const arma::uword block_rows)
{
    if (stack.n_rows <= block_rows)
    {
        return stack;
    }

    arma::mat out(stack.n_rows, stack.n_cols);
    out.head_rows(stack.n_rows - block_rows) = stack.tail_rows(stack.n_rows - block_rows);
    out.tail_rows(block_rows) = stack.tail_rows(block_rows);
    return out;
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
/// Every one of the four per-period arguments may be either one block that holds
/// in every period or a stack of one block per period, in this model's indexing:
/// block t belongs to period t. That is what serves all four dynamic factor
/// models from one function rather than from four copies of a shift convention.
///
///   - `lambda`, the measurement design: M x N, or (M tt) x N from
///     fill_stacked_loadings(). Constant is the case chan_jeliazkov_2009
///     assembles once rather than tt times, which with many observed series is
///     the difference between O(K M^2) and O(tt K M^2) on the assembly -- so a
///     drifting Lambda costs that, and there is no version of such a model that
///     does not.
///   - `u_sigma`, the measurement covariance: K x K, or (K tt) x K from
///     fill_stacked_diagonal().
///   - `v_sigma`, the covariance of the transition's innovation: N x N or
///     (N tt) x N.
///   - `a_mat`, the transition: N x Np for an order of p >= 1, or (N tt) x Np
///     from fill_stacked_transition(). A model whose factors have no dynamics
///     passes an N x N block of zeros: f_t = 0 f_{t-1} + v_t *is* the serially
///     independent -- if still heteroskedastic -- factor it then has, and its
///     prior on the first state is V, so the two cases are one code path rather
///     than a block-diagonal special case that would draw the same thing.
///
/// The last two describe the transition, so where they are stacked they go in
/// shifted; see shifted_by_one_period(). `u_sigma` and `lambda` describe the
/// measurement, whose indexing already matches, and go in as they are.
///
/// The prior over the first p states takes both transition arguments *unshifted*.
/// Those columns are the truncated transitions run from nothing, so what they
/// need is A_0 .. A_{p-1} and V_0 .. V_{p-1}, which is where the stacks start.
///
/// The last of the tt + 1 columns chan_jeliazkov_2009 returns is dropped. It is
/// the transition applied once past the end of the sample and no observation
/// informs it, so integrating it out leaves the rest exactly as it was.
inline arma::mat draw_factor_path(const arma::mat &x_t, const arma::mat &lambda,
                                  const arma::mat &u_sigma, const arma::mat &v_sigma,
                                  const arma::mat &a_mat, const int n, const int p_state)
{
    const arma::uword tt = x_t.n_cols;
    const arma::uword nn = static_cast<arma::uword>(n);
    const arma::uword prior_rows = static_cast<arma::uword>(p_state) * nn;
    const arma::vec a_init(prior_rows, arma::fill::zeros);

    const bool v_is_stacked = v_sigma.n_rows != nn;
    const bool a_is_stacked = a_mat.n_rows != nn;

    const arma::mat v_prior_blocks =
        v_is_stacked ? arma::mat(v_sigma.head_rows(prior_rows)) : v_sigma;
    const arma::mat a_prior_blocks =
        a_is_stacked ? arma::mat(a_mat.head_rows(prior_rows)) : a_mat;

    const arma::mat p_init =
        initial_state_covariance(a_prior_blocks, v_prior_blocks, n, p_state);

    const arma::mat v_transitions = v_is_stacked ? shifted_by_one_period(v_sigma, nn) : v_sigma;
    const arma::mat a_transitions = a_is_stacked ? shifted_by_one_period(a_mat, nn) : a_mat;

    return chan_jeliazkov_2009(x_t, lambda, u_sigma, v_transitions, a_transitions, a_init, p_init)
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

/// An M x N loading matrix per period, carrying nothing but the identification:
/// (M tt) x N, block t holding what identified_loadings() returns.
///
/// Built once before the chain and never rebuilt. The identifying block is the
/// half of Lambda that does not drift -- it is not drawn in any period -- so the
/// only cells fill_stacked_loadings() ever writes are the free ones, and the
/// ones and zeros put here survive the whole chain.
inline arma::mat stacked_identified_loadings(const int k, const int n, const int tt)
{
    arma::mat stack(static_cast<arma::uword>(k) * tt, n, arma::fill::zeros);
    for (int t = 0; t < tt; t++)
    {
        stack.rows(t * k, t * k + n - 1).diag().ones();
    }
    return stack;
}

/// Writes a path of free loadings into the stack chan_jeliazkov_2009 reads as a
/// measurement matrix per period: block t is Lambda_t, M x N.
///
/// `lambda` is n_lambda x tt, one period per column, each column in the row by
/// row order fill_lambda() consumes -- so this is that function once per period,
/// scattering into a submatrix instead of into a matrix.
///
/// `stack` must be (M tt) x N and must already carry the identifying block in
/// every period; stacked_identified_loadings() builds it. Nothing here writes a
/// fixed cell.
inline void fill_stacked_loadings(arma::mat &stack, const arma::mat &lambda, const int k,
                                  const int n)
{
    const int tt = static_cast<int>(lambda.n_cols);
    for (int t = 0; t < tt; t++)
    {
        int pos = 0;
        for (int i = 1; i < k; i++)
        {
            const int width = lambda_row_width(i, n);
            stack.submat(t * k + i, 0, t * k + i, width - 1) =
                arma::trans(lambda.submat(pos, t, pos + width - 1, t));
            pos += width;
        }
    }
}

/// Writes a path of transition coefficients into the stack chan_jeliazkov_2009
/// reads as a transition per period: block t is [A_1 .. A_p] of period t,
/// N x Np.
///
/// `a` is (N^2 p) x tt, one period per column, each column vec([A_1 .. A_p]) --
/// the same object DfmNormalGamma reshapes once, reshaped once per period.
///
/// Unshifted: block t is period t's own transition, which is what
/// initial_state_covariance() wants and what draw_factor_path() shifts on
/// the way into the band sampler. Shifting here instead would put the shift
/// where the residuals and the state variance also read the stack, and quietly
/// lag the model by a period in two blocks out of three.
inline void fill_stacked_transition(arma::mat &stack, const arma::mat &a, const int n,
                                    const int p)
{
    const int tt = static_cast<int>(a.n_cols);
    for (int t = 0; t < tt; t++)
    {
        stack.rows(t * n, (t + 1) * n - 1) = arma::reshape(a.col(t), n, n * p);
    }
}

/// The SUR design of the factor transition, one block per period: (N tt) x
/// (N^2 p), block t holding kron(x_t', I_N) for the lagged factors x_t of period
/// t, so that f_t = Z_t vec([A_1 .. A_p]_t).
///
/// The Kronecker product is scattered rather than formed. `a` is vec of an
/// N x Np matrix, column-major, so element (i, c) of it sits at c N + i, which
/// is the position row i of block t has to carry x_t(c) in. Every other cell is
/// structurally zero: `z_a` is expected zeroed once before the chain and the
/// cells this leaves alone stay zero for its whole life.
///
/// `x_a` is the (N p) x tt lagged factor matrix fill_lagged_factors() writes,
/// zero before the sample -- so the first p periods carry the truncated
/// transitions rather than a special case.
inline void fill_transition_design(arma::mat &z_a, const arma::mat &x_a, const int n)
{
    const int tt = static_cast<int>(x_a.n_cols);
    const int np = static_cast<int>(x_a.n_rows);
    for (int t = 0; t < tt; t++)
    {
        for (int c = 0; c < np; c++)
        {
            const double value = x_a(c, t);
            for (int i = 0; i < n; i++)
            {
                z_a(t * n + i, c * n + i) = value;
            }
        }
    }
}

/// The transition residuals of a model whose transition moves with time,
/// f_t - sum_j A_{j,t} f_{t-j}, N x tt.
///
/// The lagged factors are taken from `x_a`, which is zero before the sample, so
/// the truncation of the first p periods needs no special case -- the same
/// convention transition_residuals() follows for a constant transition, and the
/// same one the prior over the first p states is derived under.
inline arma::mat transition_residuals_tvp(const arma::mat &factors, const arma::mat &a_stack,
                                          const arma::mat &x_a, const int n)
{
    const int tt = static_cast<int>(factors.n_cols);
    arma::mat v = factors;
    for (int t = 0; t < tt; t++)
    {
        v.col(t) -= a_stack.rows(t * n, (t + 1) * n - 1) * x_a.col(t);
    }
    return v;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_DFM_SUPPORT_H
