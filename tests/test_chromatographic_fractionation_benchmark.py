"""Test regression and physics benchmark for Chromatographic Isotope Fractionation."""

from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest


def _check_weighted_errors(
    sim_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    weights: np.ndarray,
    col_pairs: list[tuple[str, str]],
    max_err: float,
    label: str,
) -> None:
    """Assert grid-weighted normalized mean error stays within max_err for all column pairs."""
    failing = {}
    for sim_col, ref_col in col_pairs:
        c_sim = sim_df[sim_col].to_numpy()
        c_ref = ref_df[ref_col].to_numpy()

        valid_mask = ~np.isnan(c_ref) & ~np.isnan(c_sim)
        if not np.any(valid_mask):
            continue

        scale = np.max(np.abs(c_ref[valid_mask]))
        if scale < 1e-12:
            continue

        normalized_pointwise_err = np.abs(c_sim[valid_mask] - c_ref[valid_mask]) / scale
        w_sub = weights[valid_mask]
        w_sub_norm = w_sub / np.sum(w_sub)
        mean_err = float(np.sum(normalized_pointwise_err * w_sub_norm))
        if mean_err > max_err:
            failing[sim_col] = mean_err

    if failing:
        report = "\n".join(
            f"  {col}: {err * 100:.3f}% error (exceeds {max_err * 100:.3f}% threshold)"
            for col, err in failing.items()
        )
        pytest.fail(f"{label} failed normalized error check:\n{report}")


def test_chromatographic_fractionation_simulation_benchmark(tmp_path: Path):
    """Run chromatographic fractionation simulation in an isolated temp directory and verify against reference & analytical solution.

    Requirements:
    1. Copies src/mantis/examples/chromatographic_fractionation to tmp_path and runs run_chromatographic_fractionation.py.
    2. Exits cleanly with returncode == 0 in under 60 seconds (typically < 2s).
    3. Compares simulated concentrations and isotope deltas with run_chromatographic_fractionation_reference.csv (<= 1%).
    4. Compares simulated concentrations and isotope deltas with analytical_ogata_banks.csv (<= 3%).
    """
    repo_root = Path(__file__).resolve().parent.parent
    src_sim_dir = repo_root / "src" / "mantis" / "examples" / "chromatographic_fractionation"
    reference_csv = src_sim_dir / "run_chromatographic_fractionation_reference.csv"
    analytical_csv = src_sim_dir / "analytical_ogata_banks.csv"

    assert src_sim_dir.is_dir(), f"Simulation directory not found: {src_sim_dir}"
    assert reference_csv.exists(), f"Reference CSV not found: {reference_csv}"
    assert analytical_csv.exists(), f"Analytical CSV not found: {analytical_csv}"

    work_dir = tmp_path / "chromatographic_fractionation"
    shutil.copytree(src_sim_dir, work_dir, ignore=shutil.ignore_patterns("__pycache__"))
    sim_script = work_dir / "run_chromatographic_fractionation.py"
    output_csv = work_dir / "run_chromatographic_fractionation.csv"
    output_csv.unlink(missing_ok=True)

    assert sim_script.exists(), f"Simulation script not found: {sim_script}"

    # 1 & 2. Run simulation with a 60-second timeout
    try:
        proc = subprocess.run(
            [sys.executable, str(sim_script)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        pytest.fail(f"Simulation timed out (> 60 seconds): {exc}")

    assert proc.returncode == 0, (
        f"Simulation exited with return code {proc.returncode}.\n"
        f"--- STDOUT ---\n{proc.stdout}\n"
        f"--- STDERR ---\n{proc.stderr}\n"
    )

    assert output_csv.exists(), f"Output CSV was not produced: {output_csv}"

    # 3. Load results, reference, and analytical benchmark
    sim_df = pd.read_csv(output_csv)
    ref_df = pd.read_csv(reference_csv)
    ana_df = pd.read_csv(analytical_csv)

    assert len(sim_df) == len(ref_df), (
        f"Grid size mismatch: sim has {len(sim_df)} points, ref has {len(ref_df)} points."
    )
    assert len(sim_df) == len(ana_df), (
        f"Grid size mismatch: sim has {len(sim_df)} points, analytical has {len(ana_df)} points."
    )

    # Calculate grid spacing weights (normalized by total domain length)
    z = ref_df["z"].to_numpy()
    dz = np.diff(z)
    dz = np.append(dz, dz[-1])
    weights = dz / (z[-1] - z[0])

    # 4. Regression check against saved reference CSV (<= 1%)
    ref_pairs = [
        (c, c)
        for c in ref_df.columns
        if (c.startswith("c_") or c.startswith("d56Fe_")) and c in sim_df.columns
    ]
    assert len(ref_pairs) > 0, "No concentration or delta columns found to compare."
    _check_weighted_errors(
        sim_df,
        ref_df,
        weights,
        ref_pairs,
        max_err=0.01,
        label="Reference regression check",
    )

    # 5. Physics check against analytical Ogata-Banks retarded advection-dispersion solution (<= 3%)
    ana_pairs = [
        ("c_Fe2_total", "Fe2_total"),
        ("c_Fe2_liq", "Fe2_liq"),
        ("c_Fe2_sorbed", "Fe2_sorbed"),
        ("d56Fe_Fe2_total", "d56Fe_Fe2_total"),
        ("d56Fe_Fe2_liq", "d56Fe_Fe2_liq"),
        ("d56Fe_Fe2_sorbed", "d56Fe_Fe2_sorbed"),
    ]
    _check_weighted_errors(
        sim_df,
        ana_df,
        weights,
        ana_pairs,
        max_err=0.03,
        label="Analytical Ogata-Banks comparison check",
    )
