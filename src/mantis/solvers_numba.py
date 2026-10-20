"""
Numba-accelerated numerical kernels for direct assembled reactive transport solvers.
"""

from __future__ import annotations

import numba as nb
import numpy as np


@nb.njit(fastmath=True)
def solve_tridiagonal_thomas(
    ab: np.ndarray,
    b: np.ndarray,
    out: np.ndarray,
    c_prime: np.ndarray,
) -> None:
    """
    O(N) In-place Thomas algorithm (TDMA) for tridiagonal systems in SciPy banded (1, 1) layout.

    Matrix layout in ab (shape 3, N):
      - row 0: superdiagonal (ab[0, 1:] contains A[i, i+1] for i = 0..N-2)
      - row 1: main diagonal  (ab[1, :] contains A[i, i] for i = 0..N-1)
      - row 2: subdiagonal   (ab[2, :-1] contains A[i+1, i] for i = 0..N-2)

    b:       RHS vector of shape (N,)
    out:     solution vector of shape (N,), updated in-place
    c_prime: work scratch array of shape (N - 1,)
    """
    n = len(b)
    denom0 = ab[1, 0]
    c_prime[0] = ab[0, 1] / denom0
    out[0] = b[0] / denom0

    # Forward elimination
    for i in range(1, n - 1):
        ai = ab[2, i - 1]
        denom = ab[1, i] - ai * c_prime[i - 1]
        c_prime[i] = ab[0, i + 1] / denom
        out[i] = (b[i] - ai * out[i - 1]) / denom

    ai_last = ab[2, n - 2]
    denom_last = ab[1, n - 1] - ai_last * c_prime[n - 2]
    out[n - 1] = (b[n - 1] - ai_last * out[n - 2]) / denom_last

    # Back substitution
    for i in range(n - 2, -1, -1):
        out[i] -= c_prime[i] * out[i + 1]


@nb.njit
def compute_species_residual_wrms(
    curr_val: np.ndarray,
    prev_val: np.ndarray,
    inner_tol: float,
    atol: float,
) -> float:
    """
    Computes WRMS norm for a single species:
        sqrt( (1/N) * sum( ( |c_curr - c_prev| / (inner_tol * |c_curr| + atol) )^2 ) )
    """
    n = len(curr_val)
    sum_sq = 0.0
    for i in range(n):
        c = curr_val[i]
        p = prev_val[i]
        if np.isnan(c) or np.isnan(p) or np.isinf(c) or np.isinf(p):
            return np.inf
        diff = abs(c - p)
        scale = inner_tol * abs(c) + atol
        if scale <= 0.0:
            return np.inf
        ratio = diff / scale
        if np.isnan(ratio) or np.isinf(ratio):
            return np.inf
        sum_sq += ratio * ratio
    return np.sqrt(sum_sq / n)


@nb.njit
def compute_species_residual_linf(
    curr_val: np.ndarray,
    prev_val: np.ndarray,
    inner_tol: float,
    atol: float,
) -> float:
    """
    Computes L-infinity norm for a single species:
        max( |c_curr - c_prev| / (inner_tol * |c_curr| + atol) )
    """
    n = len(curr_val)
    max_ratio = 0.0
    for i in range(n):
        c = curr_val[i]
        p = prev_val[i]
        if np.isnan(c) or np.isnan(p) or np.isinf(c) or np.isinf(p):
            return np.inf
        diff = abs(c - p)
        scale = inner_tol * abs(c) + atol
        if scale <= 0.0:
            return np.inf
        ratio = diff / scale
        if np.isnan(ratio) or np.isinf(ratio):
            return np.inf
        if ratio > max_ratio:
            max_ratio = ratio
    return max_ratio


@nb.njit
def compute_fractional_trust_region_alpha(
    u_curr: np.ndarray,
    delta_u: np.ndarray,
    scale_vec: np.ndarray,
    atol_vec: np.ndarray,
    max_fractional_step: float = 0.5,
) -> float:
    """
    Computes a fractional step limit alpha in (0, 1] to prevent excessive
    relative state changes in a single Newton iteration:
      - Depletion (delta_u < 0): bounded by max_fractional_step * (max(u, 0) + atol)
      - Growth (delta_u > 0):
          * if established (u >= atol): bounded by max_fractional_step * (u + atol)
          * if empty cell (u < atol): bounded by max_fractional_step * scale_vec
    """
    N, S = u_curr.shape
    alpha = 1.0
    for s in range(S):
        atol_s = atol_vec[s]
        sc_s = scale_vec[s]
        for i in range(N):
            du = delta_u[i, s]
            if abs(du) > 1e-30:
                c_val = u_curr[i, s]
                if du < 0.0:
                    allowed = max_fractional_step * (max(c_val, 0.0) + atol_s)
                    ratio = allowed / (-du)
                else:
                    base_val = max(c_val, sc_s) if c_val < atol_s else (c_val + atol_s)
                    allowed = max_fractional_step * base_val
                    ratio = allowed / du
                if ratio < alpha:
                    alpha = ratio
    return max(alpha, 1e-4)


@nb.njit
def compute_weighted_residual_norm_wrms(
    R: np.ndarray,
    inv_dt_diag: np.ndarray,
    u_curr: np.ndarray,
    atol_vec: np.ndarray,
    inner_tol: float,
) -> float:
    """
    Computes component-weighted RMS norm of residual R:
        sqrt( (1 / (N * S)) * sum_{i, s} ( R_{i, s} / (inv_dt_diag_{i, s} * (inner_tol * |u_{i, s}| + atol_s)) )^2 )
    """
    N, S = R.shape
    sum_sq = 0.0
    for s in range(S):
        atol_s = atol_vec[s]
        for i in range(N):
            r_val = R[i, s]
            if np.isnan(r_val) or np.isinf(r_val):
                return np.inf
            scale = inv_dt_diag[i, s] * (inner_tol * abs(u_curr[i, s]) + atol_s)
            if scale <= 1e-30:
                scale = 1e-30
            ratio = r_val / scale
            sum_sq += ratio * ratio
    return np.sqrt(sum_sq / (N * S))



@nb.njit(fastmath=True)
def solve_block_tridiagonal_thomas(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    out: np.ndarray,
    C_prime: np.ndarray,
    D_prime: np.ndarray,
) -> None:
    """
    O(N * S^3) In-place Block Tridiagonal Matrix Algorithm (Block TDMA) for 1D monolithic systems.

    Solves the block-tridiagonal linear system:
        A_i * U_{i-1} + B_i * U_i + C_i * U_{i+1} = D_i   for i = 0 .. N-1

    where:
      - N is the number of 1D spatial cells
      - S is the number of coupled chemical species
      - A has shape (N, S): diagonal elements of subdiagonal spatial block A_i (coupling to cell i-1)
      - B has shape (N, S, S): dense diagonal block for cell i (transport + chemical Jacobian)
      - C has shape (N, S): diagonal elements of superdiagonal spatial block C_i (coupling to cell i+1)
      - D has shape (N, S): RHS / residual vector at each cell
      - out has shape (N, S): solution vector, updated in-place
      - C_prime has shape (N, S, S): scratch array for forward block elimination
      - D_prime has shape (N, S): scratch array for forward RHS elimination
    """
    N, S = D.shape
    M = np.empty((S, S), dtype=np.float64)
    piv = np.empty(S, dtype=np.int64)
    y = np.empty(S, dtype=np.float64)
    rhs = np.empty(S, dtype=np.float64)

    # Forward elimination sweep
    for i in range(N):
        # 1. Form M_i = B_i - A_i * C'_{i-1}
        for r in range(S):
            for c in range(S):
                M[r, c] = B[i, r, c]
        if i > 0:
            for r in range(S):
                ar = A[i, r]
                for c in range(S):
                    M[r, c] -= ar * C_prime[i - 1, r, c]

        # 2. Form rhs = D_i - A_i * D'_{i-1}
        for r in range(S):
            rhs[r] = D[i, r]
        if i > 0:
            for r in range(S):
                rhs[r] -= A[i, r] * D_prime[i - 1, r]

        # 3. LU factorization of M with partial pivoting
        for k in range(S):
            piv[k] = k
        for k in range(S):
            max_r = k
            max_v = abs(M[k, k])
            for r in range(k + 1, S):
                v = abs(M[r, k])
                if v > max_v:
                    max_v = v
                    max_r = r
            if max_r != k:
                tmp_p = piv[k]
                piv[k] = piv[max_r]
                piv[max_r] = tmp_p
                for col in range(S):
                    tmp = M[k, col]
                    M[k, col] = M[max_r, col]
                    M[max_r, col] = tmp
            diag = M[k, k]
            if abs(diag) > 1e-30:
                for r in range(k + 1, S):
                    M[r, k] /= diag
                    f = M[r, k]
                    for col in range(k + 1, S):
                        M[r, col] -= f * M[k, col]

        # 4. Solve M * D'_i = rhs
        # Forward substitution Ly = P * rhs
        for r in range(S):
            val = rhs[piv[r]]
            for col in range(r):
                val -= M[r, col] * y[col]
            y[r] = val
        # Back substitution Ux = y
        for r in range(S - 1, -1, -1):
            val = y[r]
            for col in range(r + 1, S):
                val -= M[r, col] * D_prime[i, col]
            D_prime[i, r] = val / (M[r, r] if abs(M[r, r]) > 1e-30 else 1e-30)

        # 5. If i < N - 1, solve M * C'_i[:, j] = C[i, j] * e_j for each column j
        if i < N - 1:
            for j in range(S):
                cj = C[i, j]
                # Forward Ly = P * (cj * e_j)
                for r in range(S):
                    val = cj if piv[r] == j else 0.0
                    for col in range(r):
                        val -= M[r, col] * y[col]
                    y[r] = val
                # Back Ux = y
                for r in range(S - 1, -1, -1):
                    val = y[r]
                    for col in range(r + 1, S):
                        val -= M[r, col] * C_prime[i, col, j]
                    C_prime[i, r, j] = val / (M[r, r] if abs(M[r, r]) > 1e-30 else 1e-30)

    # Back substitution sweep
    for r in range(S):
        out[N - 1, r] = D_prime[N - 1, r]
    for i in range(N - 2, -1, -1):
        for r in range(S):
            val = D_prime[i, r]
            for col in range(S):
                val -= C_prime[i, r, col] * out[i + 1, col]
            out[i, r] = val

