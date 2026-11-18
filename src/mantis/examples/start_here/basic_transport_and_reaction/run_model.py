#!/usr/bin/env python
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
