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
"""MANTIS: A Monolithic Assembly of Nonlinear Transport Interactions Solver."""

import sys

import mantis.plot_data_new as _plot_data_new
from mantis.diff_lib import get_flux, make_grid_from_sections
from mantis.errors import (
    FilePathError,
    InputFileError,
    MantisError,
    MissingDefinitionError,
    ReactionParseError,
    SolverError,
    UnitParseError,
    ValidationError,
    install_excepthook,
)
from mantis.model import DiageneticModel, DiageneticModelResult
from mantis.plot_data_new import (
    isotopes,
    liquid_species,
    reaction_rates,
    solid_species,
)
from mantis.report_generator import generate_simulation_report

install_excepthook()

# Alias plot_data_new in sys.modules so 'import plot_data_new' and
# 'from plot_data_new import ...' work seamlessly even if mantis is installed as a package
sys.modules.setdefault("plot_data_new", _plot_data_new)

__all__ = [
    "DiageneticModel",
    "DiageneticModelResult",
    "generate_simulation_report",
    "get_flux",
    "make_grid_from_sections",
    "liquid_species",
    "solid_species",
    "reaction_rates",
    "isotopes",
    "MantisError",
    "InputFileError",
    "MissingDefinitionError",
    "ValidationError",
    "ReactionParseError",
    "UnitParseError",
    "FilePathError",
    "SolverError",
]

