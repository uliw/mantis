#!/usr/bin/env python
"""run_luff_model.py

Simulation of Hydrate Ridge cold vent sediments (Luff and Wallmann, 2003)
using mantis-rtm.
"""

from mantis import DiageneticModel

if __name__ == "__main__":
    model = DiageneticModel(
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        layout="plot_layout.py",
    )

    model.configure(
        grid_setup=[
            ["0 cm", "0.5 cm", "0.1 mm", "0.1 mm"],
            ["0.5 cm", "5.0 cm", "0.1 mm", "0.5 mm"],
            ["5.0 cm", "15.0 cm", "0.5 mm", "1.0 mm"],
        ],
        t_end="5 yr",
        dt_max="1 yr",
        temp=[4.2, 4.2],
        w="0.0275 cm/yr",
        advection="10 cm/yr",
        phi=0.74,
        use_numerical_jacobian=True,
        bioturbation={
            "BT0": "0.01 cm^2/year",
            "depth": "3.0 cm",
            "attenuation": "1.0 cm",
        },
    )

    print("Running Luff (2003) model simulation...")
    result = model.run()
    print("Simulation finished successfully.")
    df, fqfn = result.save_csv()
    result.save_state()
    # result.plot_summary() looks ugly
    # run plot_comparison.py instead
