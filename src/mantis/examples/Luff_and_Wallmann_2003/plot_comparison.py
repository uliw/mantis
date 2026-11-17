#!/usr/bin/env python
"""plot_comparison.py

Plot model simulation results against measured data from Station SO143/55-2
(Luff and Wallmann 2003, GCA 67: 3403-3421).

Only species/properties that exist both in the model and in the measured dataset
(Appendix, Station SO143/55-2) are plotted:
  1. Sulfate (SO4)
  2. Dissolved Calcium (Ca2+)
  3. Total Alkalinity (TA)
  4. Total Sulfide (H2S / TS2)
  5. Phosphate (TPO4)
  6. Ammonia (TNH4)
  7. Carbonate solid content (CaCO3 wt%)
  8. Porewater pH (in-situ simulated vs. ex-situ shipboard measured)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mantis.diff_lib import make_grid_from_sections, compute_carbonate_speciation

DIR = Path(__file__).resolve().parent

# 1. Grid geometry
grid_setup = [
    ["0 cm", "0.5 cm", "0.1 mm", "0.1 mm"],
    ["0.5 cm", "5.0 cm", "0.1 mm", "0.5 mm"],
    ["5.0 cm", "15.0 cm", "0.5 mm", "1.0 mm"],
]
mesh, z_m = make_grid_from_sections(grid_setup)
z_cm = np.array(mesh.cellCenters[0]) * 100.0

# 2. Load model output and measured data
state = np.load(DIR / "run_luff_model_state.npz")
meas = pd.read_csv(DIR / "data_so143_55_2.csv")

# 3. Compute in-situ speciation and pH across all cells
class Container:
    pass

class MP:
    is_carbonate = True
    K1_carb = 8.446e-7
    K2_carb = 4.533e-10
    species_defs = {
        "TS2": {"phase": "liquid", "include": True},
        "HS": {
            "phase": "liquid",
            "include": False,
            "equilibrium_of": "TS2",
            "fraction": "Ka1 / (Hplus + Ka1)",
            "constants": {"pKa1": 7.0},
            "include_in_TA": 1.0,
        },
    }

c = Container()
c.TA = state["TA"]
c.DIC = state["DIC"]
c.TS2 = state["TS2"]
mp = MP()
compute_carbonate_speciation(c, mp)
sim_ph = c.pH

# 4. Total CaCO3 in wt%
# Solid phase in mantis is in mol/m^3_solid.
# mw_CaCO3 = 100.086 g/mol, dry bulk density rho_s = 2.5 g/cm^3
# wt% = C_solid [mol/m^3] * 100.086 / (2.5 * 10000)
sim_caco3_wt = (state["Calcite"] + state["Aragonite"]) * 100.086 / 25000.0

# 5. Create 2x4 Multi-panel Plot
fig, axs = plt.subplots(2, 4, figsize=(16, 9), sharey=True)

# Panel 1: Sulfate (SO4)
axs[0, 0].plot(state["SO4"], z_cm, "b-", lw=2.2, label="Model")
axs[0, 0].plot(meas["SO4_mM"], meas["depth_cm"], "b*", markersize=9, label="Measured")
axs[0, 0].set_xlabel("SO4 [mM]")
axs[0, 0].set_ylabel("Depth [cm]")
axs[0, 0].set_ylim(15, 0)
axs[0, 0].set_xlim(-1, 30)
axs[0, 0].legend(loc="upper right")
axs[0, 0].grid(True, linestyle="--", alpha=0.5)
axs[0, 0].set_title("Sulfate (SO4)")

# Panel 2: Dissolved Calcium (Ca2+)
axs[0, 1].plot(state["Ca2"], z_cm, "g-", lw=2.2, label="Model")
axs[0, 1].plot(meas["Ca_mM"], meas["depth_cm"], "g*", markersize=9, label="Measured")
axs[0, 1].set_xlabel("Ca2+ [mM]")
axs[0, 1].set_xlim(2, 11)
axs[0, 1].legend(loc="upper right")
axs[0, 1].grid(True, linestyle="--", alpha=0.5)
axs[0, 1].set_title("Dissolved Calcium (Ca2+)")

# Panel 3: Total Alkalinity (TA)
axs[0, 2].plot(state["TA"], z_cm, "m-", lw=2.2, label="Model")
axs[0, 2].plot(meas["TA_mM"], meas["depth_cm"], "m*", markersize=9, label="Measured")
axs[0, 2].set_xlabel("TA [mM]")
axs[0, 2].set_xlim(0, 35)
axs[0, 2].legend(loc="lower right")
axs[0, 2].grid(True, linestyle="--", alpha=0.5)
axs[0, 2].set_title("Total Alkalinity (TA)")

# Panel 4: Total Sulfide (H2S / TS2)
axs[0, 3].plot(state["TS2"], z_cm, color="brown", lw=2.2, label="Model")
axs[0, 3].plot(meas["H2S_mM"], meas["depth_cm"], "*", color="brown", markersize=9, label="Measured")
axs[0, 3].set_xlabel("H2S [mM]")
axs[0, 3].set_xlim(0, 25)
axs[0, 3].legend(loc="lower right")
axs[0, 3].grid(True, linestyle="--", alpha=0.5)
axs[0, 3].set_title("Total Sulfide (TH2S)")

# Panel 5: Phosphate (TPO4)
axs[1, 0].plot(state["TPO4"] * 1000.0, z_cm, color="darkorange", lw=2.2, label="Model")
axs[1, 0].plot(meas["PO4_uM"], meas["depth_cm"], "*", color="darkorange", markersize=9, label="Measured")
axs[1, 0].set_xlabel("PO4 [uM]")
axs[1, 0].set_ylabel("Depth [cm]")
axs[1, 0].set_xlim(0, 60)
axs[1, 0].legend(loc="upper right")
axs[1, 0].grid(True, linestyle="--", alpha=0.5)
axs[1, 0].set_title("Phosphate (TPO4)")

# Panel 6: Ammonia (TNH4)
axs[1, 1].plot(state["TNH4"] * 1000.0, z_cm, color="teal", lw=2.2, label="Model")
axs[1, 1].plot(meas["NH4_uM"], meas["depth_cm"], "*", color="teal", markersize=9, label="Measured")
axs[1, 1].set_xlabel("NH4 [uM]")
axs[1, 1].set_xlim(80, 240)
axs[1, 1].legend(loc="lower right")
axs[1, 1].grid(True, linestyle="--", alpha=0.5)
axs[1, 1].set_title("Ammonia (TNH4)")

# Panel 7: Solid CaCO3 (wt%)
axs[1, 2].plot(sim_caco3_wt, z_cm, color="purple", lw=2.2, label="Model (current)")
axs[1, 2].plot(meas["CaCO3_wt_pct"], meas["depth_cm"], "*", color="purple", markersize=9, label="Measured")
axs[1, 2].set_xlabel("CaCO3 [wt%]")
# axs[1, 2].set_xlim(0, 30)
axs[1, 2].legend(loc="upper right")
axs[1, 2].grid(True, linestyle="--", alpha=0.5)
axs[1, 2].set_title("Carbonate Solid (CaCO3)")

# Panel 8: Porewater pH
axs[1, 3].plot(sim_ph, z_cm, color="crimson", lw=2.2, label="Model (in-situ)")
axs[1, 3].plot(meas["pH"], meas["depth_cm"], "*", color="crimson", markersize=9, label="Measured (degassed)")
axs[1, 3].set_xlabel("pH")
axs[1, 3].set_xlim(6.0, 9.0)
axs[1, 3].legend(loc="lower right")
axs[1, 3].grid(True, linestyle="--", alpha=0.5)
axs[1, 3].set_title("Porewater pH")

plt.tight_layout()
out_pdf = DIR / "luff_measured_vs_modeled.pdf"
plt.savefig(out_fig, dpi=200)
plt.savefig(out_pdf)
print(f"Comparison plot saved to {out_fig} and {out_pdf}")
