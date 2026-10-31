#!/usr/bin/env python
"""run_simplefe.py

Replication of the SimpleFe Diagenetic Benchmark Model (Du 2023).
==================================================================

Scientific Overview:
--------------------
This model reproduces the benchmark reactive-transport simulation "SimpleFe" described
by J. Du (2023, Geosci. Model Dev., 16, 5865-5894; https://doi.org/10.5194/gmd-16-5865-2023).
The SimpleFe benchmark evaluates the performance of diagenetic reactive-transport solvers
in simulating coupled organic carbon degradation, iron reduction, sulfate reduction,
iron monosulfide (FeS) precipitation, and porewater pH dynamics without the buffering of
solid calcium carbonate (CaCO3).

Model Architecture in Mantis:
-----------------------------
Mantis implements the SimpleFe model through declarative configuration modules:
  1. `species.py`: Chemical inventory, physical phase classification (solid, dissolved, mixed),
     molecular diffusion coefficients, and boundary conditions at the sediment-water interface.
  2. `simplefe_reactions.py`: Geochemical reaction network (DIR, OSR, sulfide oxidation,
     and FeS precipitation) with symbolic Monod/inhibition limiters.
  3. `kinetic_constants.py`: Thermodynamic equilibrium and kinetic rate constants with automatic
     unit conversion to SI base units.
  4. `plot_layout.py`: Multi-panel visualization layout comparing Mantis output directly
     against Julia SedTrace reference data (`SimpleFe_reference.csv`).

Key Physics & Geochemical Formulations:
---------------------------------------
  • Domain & Geometry: 50 cm sediment column, discretized with an adaptive geometric mesh
    refining down to 0.2 mm near the active redox transition zone (0.1 - 2.0 cm).
  • Environmental Conditions: Water depth = 500 m (~50 bar), Salinity = 35 PSU, Temp = 5 °C.
  • Diffusive Boundary Layer (DBL): Stagnant benthic boundary layer delta = 0.05 cm (0.5 mm)
    overlying the sediment-water interface, matching Du (2023) Table 6. Transport across the DBL
    is solved as an explicit series resistance: 1 / D_top = (z0 / (phi0 * D_sed)) + (delta / D_free).
    This controls the retention vs. ventilation of dissolved iron (Fe2+) and sulfide (H2S),
    accurately setting the total authigenic FeS precipitation inventory and burial flux.
  • Compaction & Porosity: Exponential porosity profile:
        phi(z) = phi_bottom + (phi_top - phi_bottom) * exp(-z / z_decay)
    with phi_top = 0.90, phi_bottom = 0.80, and decay depth = 2.0 cm (half-drop depth = 1.386 cm).
  • Burial Advection: Solid sediment burial velocity is coupled to porosity compaction to conserve
    solid mass: w_s(z) = w_inf * (1 - phi_inf) / (1 - phi(z)), with w_inf = 0.1404 cm/yr.
  • Bioturbation: Macrofaunal particle mixing follows an exponential decay with depth:
        Db(z) = Db0 * exp(-z / x_atten)
    with Db0 = 10 cm^2/yr and e-folding attenuation depth x_atten = 3.0 cm.
  • Dynamic Porewater pH: Modeled through total dissolved inorganic carbon (DIC) and total
    alkalinity (TA) transport. Carbonate speciation (Lueker et al. 2000, Millero 1995) and
    bisulfide alkalinity (Ka1 = 10^-7.0) are solved self-consistently at each grid cell.
  • Iron Speciation & Adsorption: Total iron (Fe2_total) includes both dissolved complexes
    (free Fe2+, FeHS+, FeCl+, FeSO4) and iron adsorbed onto POC (Fe_ads = K_ads * POC * Fe_dis).
    Free Fe2+ (Fe2_liq) is dynamically partitioned and drives reversible FeS precipitation.

Usage:
------
    python run_simplefe.py
    python run_simplefe.py --video "10 yr"
    python run_simplefe.py --solver-tol 1e-4
"""

from pathlib import Path
from mantis import DiageneticModel


def main():
    base_dir = Path(__file__).resolve().parent

    # 1. Initialize DiageneticModel using declarative configuration files
    model = DiageneticModel(
        base_dir=base_dir,
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        equations="equations.py",
        layout="plot_layout.py",
    )

    # 2. Configure model parameters matching Du (2023) Table 6
    model.configure(
        # video_dt = "1 year",
        grid_setup=[ # start, stop, start_dz, stop_dz
            ["0 cm", "5 cm", "50 um", "2 mm"],
            ["5 cm", "50 cm", "2 mm", "2.4 cm"],
        ],
        t_end="2 kyr",
        dt_max="10 yr",
        temp=[5.0, 5.0],
        w="0.1404 cm/yr",      # Deep burial velocity: w_inf = Fsed / rho_s / (1 - phi_bottom)

        # Diffusive Boundary Layer (DBL): in Du (2023) Table 1 & Table 6, DBL delta = 0.05 cm
        # is only specified for SO4 (which has Robin BC with DBL in SedTrace).
        # TCO2, TH2S, and TFe have Dirichlet BCs at the sediment-water interface without DBL.
        diffusive_boundary_layer=None,
        phi={
            "top": 0.9,
            "bottom": 0.8,
            "decay": "2 cm",   # Exponential decay: phi(z) = phi_bottom + (phi_top - phi_bottom) * exp(-z / decay)
            "half_depth": "1.3863 cm",  # Equivalent half-drop depth: z_half = decay * ln(2)
        },
        isotopes=False,

        # Bioturbation: Dbt(z) = 10 * exp(-z / 3 cm) cm^2/yr
        bioturbation={
            "BT0": "10 cm^2/year",
            "profile": "exponential",
            "attenuation": "3 cm",
        },
    )

    # Parse CLI flags (e.g. --video, --solver, --solver-tol)
    model.parse_cli_args()

    # 3. Execute simulation
    print("\n[SimpleFe] Running diagenesis simulation...")
    result = model.run()
    print("[SimpleFe] Simulation completed successfully.")

    # 4. Save results & render summary plot
    df, fqfn = result.save_csv()
    result.save_state()
    return result


if __name__ == "__main__":
    main()
