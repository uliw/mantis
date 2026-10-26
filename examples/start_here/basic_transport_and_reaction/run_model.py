#!/usr/bin/env python
"""A generic model file

Usage:
------
    python run_velde_model.py                       # Run standard simulation
    python run_velde_model.py --help                # Show CLI options
    python run_velde_model.py --video "10 year"     # Enable live video recording
"""

from mantis import DiageneticModel
from mantis.diff_lib import get_total_delta

if __name__ == "__main__":
    model = DiageneticModel(
        reactions="reactions.py",
        equations="equations.py",
        constants="kinetic_constants.py",
        layout="plot_layout.py",
    )

    # 2. Configure environmental and run parameters
    model.configure(
        grid_setup=[  # start, stop, start_dz, stop_dz
            ["0 cm", "100 cm", "2 mm", "2 mm"],
        ],
        t_end="2 kyr", # max run time 
        dt_max="1 year", # max dt 
        temp=[10.0, 10.1], # temperature field 
        w="0.2 cm/yr", # sedimentation rate 
        phi=0.8, # porosity 
    )

    # Allow CLI argument overrides (e.g. --video, --solver-tol)
    model.parse_cli_args()

    # 3. Run simulation
    result = model.run()

    # 4. Save results & render summary plot
    df, fqfn = result.save_csv()
    # result.save_state() optional 

    result.plot_summary()
