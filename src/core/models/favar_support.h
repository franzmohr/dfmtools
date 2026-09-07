// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_FAVAR_SUPPORT_H
#define BAYESTS_CORE_MODELS_FAVAR_SUPPORT_H

#include "bayests/arma.h"
#include "core/models/dfm_support.h"

#include <algorithm>

namespace bayests::core
{

/// What a factor augmented VAR adds to a dynamic factor model, and it is one
/// thing: half of the state is data.
///
/// Everything else is `dfm_support.h`'s -- the state path is drawn by
/// `draw_conditional_factor_path()`, the transition residuals and the diagonal
/// precision draw are the same functions, and the zero-before-the-sample
/// convention is the same convention. What lives here is the identification,
/// which is the only place the two models' arithmetic actually differs.
///
/// Throughout, `n` is the number of unobserved factors, `n_obs` the number of
/// observed ones and `n_state = n + n_obs` the width of the state vector. The
/// state is ordered factors first: s_t = (f_t', y_t')'. That order is not a
/// convention that could have gone the other way -- `chan_jeliazkov_2009's`
/// conditional entry point holds the *trailing* elements of every column fixed,
/// so the observed block has to be last for the conditioning to be a sub-block
/// extraction rather than an index map.

/// Free elements in row `i` of a k x n_state loading matrix: none while the row
/// is inside the identifying block, the whole width n_state after it.
///
/// The identifying block is the leading n x n of the factor columns, and it is
/// the *identity* -- not the unit lower triangle a dynamic factor model uses.
/// The observed columns of those same rows are zero beside it. So the first n
/// series of the panel are the factors plus idiosyncratic noise, exactly, and
/// they carry no free loading of any kind.
///
/// Why not the DFM's rule. A rotation F -> C F leaves the measurement unchanged
/// if Lambda_f absorbs it, so something has to rule C out. A DFM does it with a
/// unit lower triangular block *and* a diagonal V: both together admit only
/// C = I, which is the uniqueness of an LDL factorisation. A FAVAR has no second
/// half to offer -- its Q is unrestricted by design, that being the model -- and
/// a unit lower triangular block on its own leaves every unit lower triangular C
/// admissible. The result is a model that runs, produces plausible numbers, and
/// has loadings free to wander along an n(n-1)/2 dimensional ridge. The identity
/// block rules out C by itself: Lambda_f C^-1 has leading block C^-1, which is
/// the identity only at C = I.
///
/// The same argument covers the observed columns. F -> F + D Y is the other
/// transformation the state admits, and it takes the leading observed block from
/// zero to -D, so zero there rules it out and the two halves together pin
/// everything.
///
/// Both cases are a *prefix* of the state, which is what lets the
/// equation-by-equation draw take `state.rows(0, width - 1)` rather than
/// branching on which columns are free -- and rows inside the block are simply
/// skipped, having no equation to draw.
inline int favar_lambda_row_width(const int i, const int n, const int n_obs)
{
    return (i < n) ? 0 : n + n_obs;
}

/// A k x n_state loading matrix carrying nothing but its identification: the
/// identity in the leading n x n block of the factor columns, zeros everywhere
/// else -- the observed columns of those first n rows included.
inline arma::mat favar_identified_loadings(const int k, const int n, const int n_obs)
{
    arma::mat lambda(k, n + n_obs, arma::fill::zeros);
    for (int i = 0; i < std::min(k, n); i++)
    {
        lambda(i, i) = 1.0;
    }
    return lambda;
}

/// Writes the free loadings into Lambda, leaving the identifying block alone.
///
/// The vector runs row by row and, within a row, left to right -- the order the
/// equation-by-equation draw consumes it in, which is what makes each row's
/// slice of the prior contiguous. From row n on that means the factor loadings
/// of the row before its observed ones.
inline void fill_favar_lambda(arma::mat &lambda, const arma::vec &free, const int n,
                              const int n_obs)
{
    const int k = static_cast<int>(lambda.n_rows);
    int pos = 0;
    for (int i = 1; i < k; i++)
    {
        const int width = favar_lambda_row_width(i, n, n_obs);
        if (width == 0)
        {
            continue;
        }
        lambda.submat(i, 0, i, width - 1) = arma::trans(free.subvec(pos, pos + width - 1));
        pos += width;
    }
}

/// The observed factors as the sampler works with them: n_obs x tt, one period
/// per column.
///
/// `train.f_obs` arrives tt x n_obs, one period per row, which is the layout the
/// file stores and the layout an R caller has. Everything downstream wants the
/// transpose, so it is taken once per run rather than once per draw.
inline arma::mat obs_factors_by_period(const TrainData &train)
{
    return arma::trans(train.f_obs);
}

/// The whole state path, n_state x tt, from the drawn factors and the observed
/// factors.
///
/// Assembled once per draw and read by four of the five Gibbs blocks: the
/// loadings regress on it, the idiosyncratic residuals are taken against it, and
/// the transition both regresses on its lags and is differenced against it. The
/// alternative -- carrying the two halves separately and joining them at each
/// use -- was not worth the four rebuilds and the two index conventions.
inline arma::mat stacked_state(const arma::mat &factors, const arma::mat &obs)
{
    return arma::join_vert(factors, obs);
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_FAVAR_SUPPORT_H
