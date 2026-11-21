#!/usr/bin/env python
"""run_chromatographic_fractionation.py

Benchmark Model: Chromatographic Isotope Fractionation of Adsorbed Iron.
========================================================================

Simulates 1D advective-dispersive transport of dissolved ferrous iron invading
an initial sediment column devoid of Fe2+, undergoing equilibrium sorption to sediment:
    Fe2_liq <-> Fe2_sorbed  (f_dis = 0.5)

With isotope equilibrium fractionation between dissolved and sorbed phases:
    Fe2_sorbed_54 <-> Fe2_liq_54   (alpha_sorbed_54/liq_54 = 0.998)

Because the heavy isotope (56Fe) partitions preferentially into the adsorbed phase,
the lighter isotope (54Fe) is enriched in the porewater and experiences less retardation
than 56Fe. As the invasion front travels through the sediment column, 54Fe outpaces 56Fe,
producing a distinct negative delta56Fe chromatographic anomaly at the leading edge.

Usage:
------
    python run_chromatographic_fractionation.py
"""

from pathlib import Path
from mantis import DiageneticModel


def main():
    base_dir = Path(__file__).resolve().parent

    # 1. Initialize DiageneticModel
    model = DiageneticModel(
        base_dir=base_dir,
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        layout="plot_layout.py",
    )

    # 2. Configure model parameters
    # 20 cm grid, uniform 1 mm spacing
    model.configure(
        grid_setup=[
            ["0 cm", "20 cm", "1 mm", "1 mm"],
        ],
        t_end="0.5 yr",
        dt_init="0.5 day",
        dt_max="2 day",
        temp=[10.0, 10.0],
        w="0 cm/yr",
        advection="-5 cm/yr",   # Downward porewater advection: 5 cm/yr (u > 0 downward)
        diffusive_boundary_layer=None,
        phi=0.8,                # Uniform constant porosity
        isotopes=True,
        delta_threshold_liquid=1e-4,
        delta_threshold_solid=1e-4,
        bioturbation=None,
        bioirrigation=None,
    )

    # Parse CLI flags (e.g. --solver, --solver-tol)
    model.parse_cli_args()

    # 3. Run simulation
    print("\n[ChromatographicFractionation] Running simulation...")
    result = model.run()
    print("[ChromatographicFractionation] Simulation completed successfully.")

    # 4. Save results & render summary plot
    df, fqfn = result.save_csv()
    result.save_state()
    try:
        result.plot_summary(outfile=str(base_dir / "run_chromatographic_fractionation.svg"))
    except Exception as e:
        print(f"[Warning] Could not generate summary plot: {e}")
    return result


if __name__ == "__main__":
    main()
