"""Test regression benchmark for Velde et al. (2016) simulation."""

from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest


def test_velde_et_al_2016_simulation_benchmark():
    """
    Run Velde et al. (2016) simulation and verify accuracy against reference results.

    Requirements:
    1. Runs examples/velde_et_al_2016/no_isotopes/run_velde_fast.py.
    2. Exits cleanly with returncode == 0 in under 60 seconds.
    3. Compares simulated concentrations with run_velde_fast_reference.csv.
    4. Asserts normalized mean error across grid spacing is <= 1% for every active species.
    """
    repo_root = Path(__file__).resolve().parent.parent
    sim_dir = repo_root / "examples" / "velde_et_al_2016" / "no_isotopes"
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

    # 4. Compare all concentration columns (c_*)
    conc_cols = [c for c in ref_df.columns if c.startswith("c_") and c in sim_df.columns]
    assert len(conc_cols) > 0, "No concentration columns found to compare."

    errors = {}
    failing_species = {}

    for col in conc_cols:
        c_ref = ref_df[col].to_numpy()
        c_sim = sim_df[col].to_numpy()
        scale = np.max(np.abs(c_ref))

        # Skip species that are identically zero (inactive or zero-background)
        if scale < 1e-15:
            continue

        normalized_pointwise_err = np.abs(c_sim - c_ref) / scale
        mean_err = np.sum(normalized_pointwise_err * weights)
        errors[col] = mean_err

        if mean_err > 0.01:
            failing_species[col] = mean_err

    if failing_species:
        report = "\n".join(
            f"  {col}: {err * 100:.3f}% error (exceeds 1.000% threshold)"
            for col, err in failing_species.items()
        )
        pytest.fail(f"Simulation accuracy check failed for species:\n{report}")
