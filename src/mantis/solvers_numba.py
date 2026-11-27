"""Mantis-RTM: A Python library for reaction-transport modeling

Copyright (C), 2026 - Ulrich G. Wortmann

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import annotations

"""Numba-accelerated numerical kernels for direct assembled reactive transport solvers."""

import ctypes
import importlib
import os
from collections.abc import Generator
from contextlib import contextmanager
from typing import Any

import numba as nb
import numpy as np
from scipy.linalg.lapack import dgbsv

BLAS_THREAD_ENV_VARS: tuple[str, ...] = (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "OMP_NUM_THREADS",
)
LARGE_SPECIES_BLAS_THRESHOLD: int = 80
MAX_AUTO_BLAS_THREADS: int = 4

_BLAS_THREAD_API_SYMBOLS: tuple[tuple[str, str], ...] = (
    ("openblas_get_num_threads", "openblas_set_num_threads"),
    ("scipy_openblas_get_num_threads64_", "scipy_openblas_set_num_threads64_"),
    ("scipy_openblas_get_num_threads", "scipy_openblas_set_num_threads"),
    ("MKL_Get_Max_Threads", "MKL_Set_Num_Threads"),
    ("bli_thread_get_num_threads", "bli_thread_set_num_threads"),
)


def _find_blas_thread_api() -> tuple[Any, Any] | None:
    """Locate runtime BLAS thread getter/setter C symbols if available."""
    candidates: list[str | None] = []
    try:
        flapack_path = getattr(importlib.import_module("scipy.linalg._flapack"), "__file__", None)
    except Exception:
        flapack_path = None
    if isinstance(flapack_path, str) and flapack_path:
        candidates.append(flapack_path)
    if os.name == "posix":
        candidates.append(None)

    for candidate in candidates:
        try:
            lib = ctypes.CDLL(candidate)
        except OSError:
            continue
        for get_name, set_name in _BLAS_THREAD_API_SYMBOLS:
            getter = getattr(lib, get_name, None)
            setter = getattr(lib, set_name, None)
            if getter is not None and setter is not None:
                getter.restype = ctypes.c_int
                setter.argtypes = [ctypes.c_int]
                return getter, setter
    return None


def normalize_blas_threads(value: Any, param_name: str = "blas_threads") -> int | None:
    """Validate and normalize a BLAS thread limit setting to ``int`` or ``None``."""
    if value is None or (isinstance(value, str) and value.strip().lower() in ("none", "auto", "")):
        return None
    if isinstance(value, bool):
        raise ValueError(f"{param_name} must be a positive integer or None")
    if isinstance(value, int):
        threads = value
    elif isinstance(value, str):
        try:
            threads = int(value.strip())
        except ValueError as exc:
            raise ValueError(f"{param_name} must be a positive integer or None") from exc
    else:
        raise ValueError(f"{param_name} must be a positive integer or None")
    if threads < 1:
        raise ValueError(f"{param_name} must be a positive integer or None")
    return threads


def recommended_blas_threads(num_species: int, requested: int | None = None) -> int | None:
    """Determine the runtime BLAS thread limit for banded LAPACK solves.

    Explicit user configuration via ``requested`` takes highest precedence,
    followed by any user-set BLAS/OpenMP environment variable (in which case
    ``None`` is returned so the environment setting is preserved). Otherwise,
    returns ``1`` for ``num_species < 80`` (avoiding thread-pool wakeup and
    oversubscription overhead) and ``min(4, os.cpu_count() or 1)`` for
    ``num_species >= 80`` (capturing panel-BLAS speedup without cross-chiplet
    contention).

    Parameters
    ----------
    num_species : int
        Number of coupled species in the monolithic system.
    requested : int or None, optional
        Explicit thread override configured by the user.

    Returns
    -------
    int or None
        Target BLAS thread limit, or ``None`` to leave the runtime default untouched.
    """
    threads = normalize_blas_threads(requested, param_name="blas_threads")
    if threads is not None:
        return threads

    for env_var in BLAS_THREAD_ENV_VARS:
        if os.environ.get(env_var, "").strip():
            return None

    if num_species < LARGE_SPECIES_BLAS_THRESHOLD:
        return 1

    return min(MAX_AUTO_BLAS_THREADS, os.cpu_count() or 1)


@contextmanager
def blas_thread_limit(num_threads: int | None) -> Generator[int | None]:
    """Temporarily limit BLAS/LAPACK worker threads within a context block.

    Parameters
    ----------
    num_threads : int or None
        Maximum BLAS thread count for the duration of the context, or ``None``
        to leave the current BLAS configuration unchanged.

    Yields
    ------
    int or None
        Active thread limit applied inside the context, or ``None`` if no limit
        was applied.
    """
    threads = normalize_blas_threads(num_threads, param_name="num_threads")
    if threads is None:
        yield None
        return

    try:
        threadpoolctl = importlib.import_module("threadpoolctl")
    except ImportError:
        threadpoolctl = None

    if threadpoolctl is not None:
        with threadpoolctl.threadpool_limits(limits=threads, user_api="blas"):
            yield threads
        return

    api = _find_blas_thread_api()
    if api is None:
        yield None
        return

    getter, setter = api
    prev_threads = int(getter())
    setter(threads)
    try:
        yield threads
    finally:
        if prev_threads > 0:
            setter(prev_threads)


@nb.njit
def compute_species_residual_wrms(
    curr_val: np.ndarray,
    prev_val: np.ndarray,
    inner_tol: float,
    atol: float,
) -> float:
    r"""Compute WRMS norm for a single species.

    .. math::

        \sqrt{\frac{1}{N} \sum \left( \frac{|c_{\text{curr}} - c_{\text{prev}}|}{\text{inner\_tol} \cdot |c_{\text{curr}}| + \text{atol}} \right)^2}
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
    r"""Compute L-infinity norm for a single species.

    .. math::

        \max \frac{|c_{\text{curr}} - c_{\text{prev}}|}{\text{inner\_tol} \cdot |c_{\text{curr}}| + \text{atol}}
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
    """Compute a fractional step limit alpha in (0, 1] to prevent excessive state changes.

    Limits relative state changes in a single Newton iteration:
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
    r"""Compute component-weighted RMS norm of residual R.

    .. math::

        \sqrt{\frac{1}{N \cdot S} \sum_{i, s} \left( \frac{R_{i, s}}{\text{inv\_dt\_diag}_{i, s} \cdot (\text{inner\_tol} \cdot |u_{i, s}| + \text{atol}_s)} \right)^2}
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


@nb.njit
def pack_block_tridiagonal_dgbsv(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    ab_f: np.ndarray,
) -> None:
    """Pack a 1D block-tridiagonal matrix into Fortran-ordered LAPACK ``dgbsv`` band storage.

    Parameters
    ----------
    A : ndarray, shape (N, S)
        Subdiagonal spatial coupling coefficients.
    B : ndarray, shape (N, S, S)
        Dense per-cell diagonal coupling blocks.
    C : ndarray, shape (N, S)
        Superdiagonal spatial coupling coefficients.
    ab_f : ndarray, shape (3 * S + 1, N * S)
        Preallocated Fortran-contiguous LAPACK band buffer (updated in-place).
    """
    N, S = A.shape
    ab_f[:, :] = 0.0
    for i in range(N):
        base = i * S
        if i > 0:
            for s in range(S):
                ab_f[3 * S, base - S + s] = A[i, s]
        if i < N - 1:
            for s in range(S):
                ab_f[S, base + S + s] = C[i, s]
        for c in range(S):
            col = base + c
            for r in range(S):
                ab_f[2 * S + r - c, col] = B[i, r, c]


@nb.njit
def pack_block_tridiagonal_dgbsv_from_jacobian(
    A: np.ndarray,
    inv_dt: np.ndarray,
    C: np.ndarray,
    vol: np.ndarray,
    J_chem: np.ndarray,
    ab_f: np.ndarray,
) -> None:
    """Pack ``diag(inv_dt) - vol * J_chem`` and spatial stencils directly into ``dgbsv`` band storage.

    Parameters
    ----------
    A : ndarray, shape (N, S)
        Subdiagonal spatial coupling coefficients.
    inv_dt : ndarray, shape (N, S)
        Diagonal accumulation + transport stencil terms.
    C : ndarray, shape (N, S)
        Superdiagonal spatial coupling coefficients.
    vol : ndarray, shape (N,)
        Cell volumes.
    J_chem : ndarray, shape (N, S, S)
        Chemical Jacobian blocks ``dR_s / dC_m``.
    ab_f : ndarray, shape (3 * S + 1, N * S)
        Preallocated Fortran-contiguous LAPACK band buffer (updated in-place).
    """
    N, S = A.shape
    ab_f[:, :] = 0.0
    for i in range(N):
        base = i * S
        neg_vi = -vol[i]
        if i > 0:
            for s in range(S):
                ab_f[3 * S, base - S + s] = A[i, s]
        if i < N - 1:
            for s in range(S):
                ab_f[S, base + S + s] = C[i, s]
        for c in range(S):
            col = base + c
            for r in range(S):
                val = neg_vi * J_chem[i, r, c]
                if r == c:
                    val += inv_dt[i, r]
                ab_f[2 * S + r - c, col] = val


def solve_dgbsv_banded(
    S: int,
    ab_f: np.ndarray,
    D: np.ndarray,
    out: np.ndarray,
) -> None:
    """Solve a packed banded block-tridiagonal system in-place via LAPACK ``dgbsv``.

    Parameters
    ----------
    S : int
        Number of coupled chemical species (lower and upper bandwidth ``kl = ku = S``).
    ab_f : ndarray, shape (3 * S + 1, N * S)
        Fortran-ordered LAPACK band matrix (overwritten in-place).
    D : ndarray, shape (N, S)
        Right-hand side vector.
    out : ndarray, shape (N, S)
        Output solution array (updated in-place).
    """
    _, _, x, info = dgbsv(S, S, ab_f, D.ravel(), overwrite_ab=1, overwrite_b=0)
    if info > 0:
        raise np.linalg.LinAlgError(f"Singular banded Jacobian matrix in dgbsv (info={info})")
    if info < 0:
        raise ValueError(f"Illegal argument {-info} passed to LAPACK dgbsv")
    out[:] = x.reshape(D.shape)


def solve_block_tridiagonal_banded(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    out: np.ndarray,
    ab_f: np.ndarray | None = None,
) -> None:
    """Solve a 1D block-tridiagonal system using LAPACK ``dgbsv`` banded LU factorization.

    Parameters
    ----------
    A : ndarray, shape (N, S)
        Subdiagonal spatial coupling coefficients.
    B : ndarray, shape (N, S, S)
        Dense per-cell diagonal coupling blocks.
    C : ndarray, shape (N, S)
        Superdiagonal spatial coupling coefficients.
    D : ndarray, shape (N, S)
        Right-hand side vector.
    out : ndarray, shape (N, S)
        Output solution array (updated in-place).
    ab_f : ndarray or None, optional
        Optional preallocated Fortran-ordered band buffer of shape ``(3 * S + 1, N * S)``.
    """
    N, S = D.shape
    if ab_f is None:
        ab_f = np.zeros((3 * S + 1, N * S), dtype=np.float64, order="F")
    pack_block_tridiagonal_dgbsv(A, B, C, ab_f)
    solve_dgbsv_banded(S, ab_f, D, out)

