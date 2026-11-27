"""
Unit tests for Numba-accelerated solvers in mantis.solvers_numba.
"""

import numpy as np
import pytest


def test_compute_species_residual_wrms_accuracy():
    """Verify compute_species_residual_wrms produces identical result to NumPy reference."""
    from mantis.solvers_numba import compute_species_residual_wrms, compute_species_residual_linf

    rng = np.random.default_rng(99)
    N = 500
    curr = rng.uniform(0.1, 10.0, N)
    prev = curr + rng.normal(0, 0.05, N)
    inner_tol = 1e-4
    atol = 1e-6

    # Reference WRMS
    diff = np.abs(curr - prev)
    scale = inner_tol * np.abs(curr) + atol
    ratio = diff / scale
    expected_wrms = float(np.sqrt(np.mean(ratio**2)))
    expected_linf = float(np.max(ratio))

    val_wrms = compute_species_residual_wrms(curr, prev, inner_tol, atol)
    val_linf = compute_species_residual_linf(curr, prev, inner_tol, atol)

    assert abs(val_wrms - expected_wrms) < 1e-12
    assert abs(val_linf - expected_linf) < 1e-12


def test_compute_fractional_trust_region_alpha():
    """Verify compute_fractional_trust_region_alpha properly bounds excessive steps."""
    from mantis.solvers_numba import compute_fractional_trust_region_alpha

    N, S = 10, 3
    u_curr = np.ones((N, S), dtype=np.float64) * 10.0
    delta_u = np.zeros((N, S), dtype=np.float64)
    scale_vec = np.array([10.0, 10.0, 10.0], dtype=np.float64)
    atol_vec = np.array([1e-6, 1e-6, 1e-6], dtype=np.float64)

    # 1. Zero change -> alpha = 1.0
    alpha = compute_fractional_trust_region_alpha(u_curr, delta_u, scale_vec, atol_vec, max_fractional_step=0.5)
    assert alpha == 1.0

    # 2. Moderate change (delta_u = -2.0, which is 20% of 10.0) -> under 50% limit -> alpha = 1.0
    delta_u[0, 0] = -2.0
    alpha = compute_fractional_trust_region_alpha(u_curr, delta_u, scale_vec, atol_vec, max_fractional_step=0.5)
    assert alpha == 1.0

    # 3. Excessive depletion (delta_u = -10.0, 100% depletion) -> limited to 0.5 * 10 / 10 = 0.5
    delta_u[0, 0] = -10.0
    alpha = compute_fractional_trust_region_alpha(u_curr, delta_u, scale_vec, atol_vec, max_fractional_step=0.5)
    assert alpha == pytest.approx(0.5, rel=1e-4)


def test_compute_weighted_residual_norm_wrms():
    """Verify compute_weighted_residual_norm_wrms matches reference NumPy calculation."""
    from mantis.solvers_numba import compute_weighted_residual_norm_wrms

    N, S = 20, 4
    rng = np.random.default_rng(123)
    R = rng.uniform(-1.0, 1.0, (N, S))
    inv_dt = rng.uniform(0.5, 2.0, (N, S))
    u_curr = rng.uniform(0.1, 10.0, (N, S))
    atol_vec = np.full(S, 1e-5)
    inner_tol = 1e-3

    # Reference
    scale = inv_dt * (inner_tol * np.abs(u_curr) + atol_vec[None, :])
    expected = float(np.sqrt(np.mean((R / scale) ** 2)))

    val = compute_weighted_residual_norm_wrms(R, inv_dt, u_curr, atol_vec, inner_tol)
    assert val == pytest.approx(expected, rel=1e-12)


def test_solve_block_tridiagonal_banded_accuracy():
    """Verify solve_block_tridiagonal_banded matches dense scipy.linalg.solve to < 1e-11."""
    from scipy.linalg import solve
    from mantis.solvers_numba import solve_block_tridiagonal_banded

    rng = np.random.default_rng(42)

    for N, S in [(3, 1), (3, 2), (5, 4), (50, 5), (100, 10), (20, 25)]:
        A = -rng.uniform(0.1, 1.0, (N, S))
        C = -rng.uniform(0.1, 1.0, (N, S))
        B = rng.uniform(-0.2, 0.2, (N, S, S))
        for i in range(N):
            for r in range(S):
                B[i, r, r] += 4.0 + abs(A[i, r]) + abs(C[i, r])
        D = rng.uniform(1.0, 5.0, (N, S))

        K = np.zeros((N * S, N * S), dtype=np.float64)
        for i in range(N):
            K[i * S : (i + 1) * S, i * S : (i + 1) * S] = B[i]
            if i > 0:
                for r in range(S):
                    K[i * S + r, (i - 1) * S + r] = A[i, r]
            if i < N - 1:
                for r in range(S):
                    K[i * S + r, (i + 1) * S + r] = C[i, r]
        expected = solve(K, D.ravel()).reshape((N, S))

        out_auto = np.zeros((N, S), dtype=np.float64)
        solve_block_tridiagonal_banded(A, B, C, D, out_auto, ab_f=None)
        assert np.max(np.abs(out_auto - expected)) < 1e-11

        ab_f = np.zeros((3 * S + 1, N * S), dtype=np.float64, order="F")
        out_prealloc = np.zeros((N, S), dtype=np.float64)
        solve_block_tridiagonal_banded(A, B, C, D, out_prealloc, ab_f=ab_f)
        assert np.max(np.abs(out_prealloc - expected)) < 1e-11


def test_pack_block_tridiagonal_dgbsv_from_jacobian_accuracy():
    """Verify fused Jacobian packer + dgbsv matches dense scipy.linalg.solve to < 1e-11."""
    from scipy.linalg import solve
    from mantis.solvers_numba import (
        pack_block_tridiagonal_dgbsv_from_jacobian,
        solve_dgbsv_banded,
    )

    rng = np.random.default_rng(101)

    for N, S in [(5, 1), (5, 2), (40, 10)]:
        A = -rng.uniform(0.1, 1.0, (N, S))
        C = -rng.uniform(0.1, 1.0, (N, S))
        inv_dt = rng.uniform(5.0, 15.0, (N, S)) + np.abs(A) + np.abs(C)
        vol = rng.uniform(0.01, 0.5, N)
        J_chem = rng.uniform(-0.5, 0.5, (N, S, S))
        D = rng.uniform(-3.0, 3.0, (N, S))

        B = -vol[:, None, None] * J_chem
        for s in range(S):
            B[:, s, s] += inv_dt[:, s]

        K = np.zeros((N * S, N * S), dtype=np.float64)
        for i in range(N):
            K[i * S : (i + 1) * S, i * S : (i + 1) * S] = B[i]
            if i > 0:
                for r in range(S):
                    K[i * S + r, (i - 1) * S + r] = A[i, r]
            if i < N - 1:
                for r in range(S):
                    K[i * S + r, (i + 1) * S + r] = C[i, r]
        expected = solve(K, D.ravel()).reshape((N, S))

        ab_f = np.zeros((3 * S + 1, N * S), dtype=np.float64, order="F")
        out_dgbsv = np.zeros((N, S), dtype=np.float64)
        pack_block_tridiagonal_dgbsv_from_jacobian(A, inv_dt, C, vol, J_chem, ab_f)
        solve_dgbsv_banded(S, ab_f, D, out_dgbsv)

        assert np.max(np.abs(out_dgbsv - expected)) < 1e-11


def test_solve_block_tridiagonal_banded_singular_raises():
    """Verify a singular matrix raises numpy.linalg.LinAlgError."""
    from mantis.solvers_numba import solve_block_tridiagonal_banded

    N, S = 4, 3
    A = np.zeros((N, S), dtype=np.float64)
    C = np.zeros((N, S), dtype=np.float64)
    B = np.broadcast_to(np.eye(S, dtype=np.float64), (N, S, S)).copy()
    B[2, 1, :] = 0.0
    D = np.ones((N, S), dtype=np.float64)
    out = np.zeros((N, S), dtype=np.float64)

    with pytest.raises(np.linalg.LinAlgError, match="Singular banded Jacobian matrix"):
        solve_block_tridiagonal_banded(A, B, C, D, out)


def test_recommended_blas_threads_policy(monkeypatch):
    """Verify automatic and explicit BLAS thread policy rules."""
    import os

    from mantis.solvers_numba import BLAS_THREAD_ENV_VARS, recommended_blas_threads

    for var in BLAS_THREAD_ENV_VARS:
        monkeypatch.delenv(var, raising=False)

    assert recommended_blas_threads(10) == 1
    assert recommended_blas_threads(79) == 1
    expected_large = min(4, os.cpu_count() or 1)
    assert recommended_blas_threads(80) == expected_large
    assert recommended_blas_threads(100) == expected_large

    # Environment variable overrides auto selection cleanly
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "2")
    assert recommended_blas_threads(10) is None
    assert recommended_blas_threads(100) is None

    # Explicit user override wins even over environment variables
    assert recommended_blas_threads(10, requested=3) == 3
    assert recommended_blas_threads(100, requested=2) == 2

    for bad in (0, -2, True):
        with pytest.raises(ValueError, match="blas_threads"):
            recommended_blas_threads(10, requested=bad)


def test_blas_thread_limit_context_manager():
    """Verify blas_thread_limit updates, solves accurately, and restores thread counts."""
    from mantis.solvers_numba import (
        _find_blas_thread_api,
        blas_thread_limit,
        solve_block_tridiagonal_banded,
    )

    with blas_thread_limit(None) as active:
        assert active is None

    for bad in (0, -1, False):
        with pytest.raises(ValueError, match="num_threads"), blas_thread_limit(bad):
            pass

    api = _find_blas_thread_api()
    before = int(api[0]()) if api is not None else None

    with blas_thread_limit(1) as active:
        if api is not None:
            assert active == 1
            assert int(api[0]()) == 1

        N, S = 8, 3
        A = np.full((N, S), -0.1, dtype=np.float64)
        C = np.full((N, S), -0.1, dtype=np.float64)
        B = np.broadcast_to(4.0 * np.eye(S, dtype=np.float64), (N, S, S)).copy()
        D = np.ones((N, S), dtype=np.float64)
        out = np.zeros((N, S), dtype=np.float64)
        solve_block_tridiagonal_banded(A, B, C, D, out)
        assert np.all(np.isfinite(out))

    if api is not None and before is not None:
        assert int(api[0]()) == before

    with pytest.raises(RuntimeError, match="boom"), blas_thread_limit(2) as active:
        if api is not None:
            assert active == 2
            assert int(api[0]()) == 2
        raise RuntimeError("boom")

    if api is not None and before is not None:
        assert int(api[0]()) == before


