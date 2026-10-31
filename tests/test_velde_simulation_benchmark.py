"""Test regression benchmark for Velde et al. (2016) simulation (no_isotopes and with_isotopes)."""

from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest


def _run_and_verify_simulation_benchmark(scenario_dir_name: str, check_isotopes: bool = False):
    """
    Run Velde et al. (2016) simulation and verify accuracy against reference results.

    Requirements:
    1. Runs examples/velde_et_al_2016/{scenario_dir_name}/run_velde_fast.py.
    2. Exits cleanly with returncode == 0 in under 60 seconds.
    3. Compares simulated concentrations (and isotope delta values if check_isotopes=True)
       with run_velde_fast_reference.csv.
    4. Asserts normalized mean error across grid spacing is <= 1% for every active species/isotope.
    """
    repo_root = Path(__file__).resolve().parent.parent
    sim_dir = repo_root / "examples" / "velde_et_al_2016" / scenario_dir_name
    sim_script = sim_dir / "run_velde_fast.py"
    output_csv = sim_dir / "run_velde_fast.csv"
    reference_csv = sim_dir / "run_velde_fast_reference.csv"

    assert sim_script.exists(), f"Simulation script not found: {sim_script}"
    assert reference_csv.exists(), f"Reference CSV not found: {reference_csv}"

    # 1 & 2. Run simulation with a 60-second timeout
    try:
        proc = subprocess.run(
            [sys.executable, str(sim_script)],
            cwd=str(sim_dir),
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

    # 3. Load results and reference
    sim_df = pd.read_csv(output_csv)
    ref_df = pd.read_csv(reference_csv)

    assert len(sim_df) == len(ref_df), (
        f"Grid size mismatch: sim has {len(sim_df)} points, ref has {len(ref_df)} points."
    )

    # Calculate grid spacing weights (normalized by total domain length)
    z = ref_df["z"].to_numpy()
    dz = np.diff(z)
    dz = np.append(dz, dz[-1])
    total_length = z[-1] - z[0]
    weights = dz / total_length

    # 4. Compare all concentration columns (c_*) and isotope delta columns (d34S_*)
    cols_to_check = [c for c in ref_df.columns if c.startswith("c_") and c in sim_df.columns]
    if check_isotopes:
        cols_to_check += [c for c in ref_df.columns if c.startswith("d34S_") and c in sim_df.columns]
    assert len(cols_to_check) > 0, "No concentration columns found to compare."

    errors = {}
    failing_species = {}

    for col in cols_to_check:
        c_ref = ref_df[col].to_numpy()
        c_sim = sim_df[col].to_numpy()

        valid_mask = ~np.isnan(c_ref) & ~np.isnan(c_sim)
        if not np.any(valid_mask):
            continue

        scale = np.max(np.abs(c_ref[valid_mask]))

        # Skip species that are identically zero (inactive or zero-background)
        if scale < 1e-15:
            continue

        normalized_pointwise_err = np.abs(c_sim[valid_mask] - c_ref[valid_mask]) / scale
        w_sub = weights[valid_mask]
        w_sub_norm = w_sub / np.sum(w_sub)
        mean_err = np.sum(normalized_pointwise_err * w_sub_norm)
        errors[col] = mean_err

        if mean_err > 0.01:
            failing_species[col] = mean_err

    if failing_species:
        report = "\n".join(
            f"  {col}: {err * 100:.3f}% error (exceeds 1.000% threshold)"
            for col, err in failing_species.items()
        )
        pytest.fail(f"Simulation accuracy check failed for {scenario_dir_name}:\n{report}")


def test_velde_et_al_2016_simulation_benchmark():
    """Run regression benchmark for no_isotopes scenario."""
    _run_and_verify_simulation_benchmark("no_isotopes", check_isotopes=False)


def test_velde_et_al_2016_no_isotopes_with_irrigation_simulation_benchmark():
    """Run regression benchmark for no_isotopes_with_irrigation scenario."""
    _run_and_verify_simulation_benchmark("no_isotopes_with_irrigation", check_isotopes=False)


def test_velde_et_al_2016_with_isotopes_simulation_benchmark():
    """Run regression benchmark for with_isotopes scenario."""
    _run_and_verify_simulation_benchmark("with_isotopes", check_isotopes=True)
