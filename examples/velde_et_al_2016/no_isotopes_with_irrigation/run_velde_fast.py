#!/usr/bin/env python
"""run_velde_model.py

10-meter Velde Sediment Diagenesis Simulation with Pyrite Precipitation and Sulfur Isotopes.

Solves the coupled reactive-transport system using the SymPy analytical Jacobian
and Numba Block TDMA solver (solve_block_tridiagonal_thomas), driven by the
declarative DiageneticModel interface.

Usage:
------
    python run_velde_model.py                       # Run standard simulation
    python run_velde_model.py --help                # Show CLI options
    python run_velde_model.py --video "10 year"     # Enable live video recording
    python run_velde_model.py --numerical-jacobian  # Use finite-difference Jacobian
"""

from mantis import DiageneticModel
from mantis.diff_lib import get_total_delta

if __name__ == "__main__":
    # 1. Initialize model (automatically detects species.py, reactions, etc. in same dir,
    # and derives experiment name from the script filename 'run_velde_model')
    model = DiageneticModel(
        reactions="reactions.py",
        equations="equations.py",
        constants="kinetic_constants_fast.py",
        layout="plot_layout.py",
    )

    # 2. Configure environmental and run parameters
    model.configure(
        # video_dt = "1 year",
        grid_setup=[
            ["0 cm", "0.1 cm", "0.1 mm", "0.05 mm"],
            ["0.1 cm", "4.0 cm", "0.05 mm", "0.05 mm"],
            ["4.0 cm", "50 cm", "0.05 mm", "1.0 cm"],
        ],
        t_end="2 kyr",
        dt_max="1 year",
        temp=[10.0, 10.1],
        w="0.2 cm/yr",
        phi=0.8,
        isotopes=False,
        bioturbation={
            "BT0": "4 cm^2/year",
            "depth": "7.63 cm",
            "attenuation": "2 cm",
        },
        bioirrigation={
            "BI0": "182.625 /year",  # a0 = 0.5 /day = 182.625 /year
            "attenuation": "3 cm",   # xirr = 3 cm from Table 1 / Fig 2b
        },
    )

    # Allow CLI argument overrides (e.g. --video, --solver-tol)
    model.parse_cli_args()

    # 3. Run simulation
    result = model.run()

    # 4. Save results & render summary plot
    df, fqfn = result.save_csv()
    result.save_state()

    total_iron = df.c_Fe2_total + df.c_Fe3 + df.c_FeS
    print(f"total_iron diff = {total_iron.max() - total_iron.min():.2e}")
    print(f"max FeS = {df.c_FeS.max():.4e} mmol/L")
    print(f"max TS2 = {df.c_TS2.max():.4e} mmol/L")
    print(f"max Fe2 = {df.c_Fe2_total.max():.4e} mmol/L")

    result.plot_summary()
