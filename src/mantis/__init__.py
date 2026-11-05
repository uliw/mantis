"""MANTIS: A Monolithic Assembly of Nonlinear Transport Interactions Solver."""

import sys
from mantis.model import DiageneticModel, DiageneticModelResult
from mantis.diff_lib import make_grid_from_sections
from mantis.plot_data_new import (
    liquid_species,
    solid_species,
    reaction_rates,
    isotopes,
)
import mantis.plot_data_new as _plot_data_new

# Alias plot_data_new in sys.modules so 'import plot_data_new' and
# 'from plot_data_new import ...' work seamlessly even if mantis is installed as a package
sys.modules.setdefault("plot_data_new", _plot_data_new)

__all__ = [
    "DiageneticModel",
    "DiageneticModelResult",
    "make_grid_from_sections",
    "liquid_species",
    "solid_species",
    "reaction_rates",
    "isotopes",
]

