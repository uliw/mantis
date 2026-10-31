"""
plot_layout.py

Declarative layout configuration for the SimpleFe benchmark model (Du 2023).
"""

# ----  optional comparison against existing data --- #
comparison_data = "SimpleFe_Du_2023.csv"

# Conversion factors:
# - rho_s = 2.6 g/cm^3 (dry solid density)
# - mw_Fe = 55.845 g/mol (molecular weight of iron)
#
# SedTrace (Du 2023) reference units:
#   • depth: cm -> meters: 1 / 100
#   • Solutes (SO4, TH2S, TFe): mmol/cm^3_pw -> mmol/L_pw (mM): * 1000.0
#   • Solids (FeOOH, FeS): mmol/cm^3_solid -> wt% Fe:
#       wt% Fe = C_solid [mmol/cm^3] * mw_Fe [g/mol] / (10 * rho_s [g/cm^3])
#              = C_solid * 55.845 / 26.0 (~2.14788)
#   • pH: direct pH units (free scale): 1.0

rho_s = 2.6
mw_Fe = 55.845
fe_to_wt = mw_Fe / (10.0 * rho_s)  # 55.845 / 26.0

# Species mapping between mantis and external SedTrace (Du 2023) reference data.
# Format: "mantis_species": [external_column_name, source_unit, target_display_unit]
#
# Reference Data Units in SimpleFe_Du_2023.csv:
#   • depth: cm (converted to m in mantis)
#   • SO4, TH2S, TFe, Fe_aq (porewater solutes): mmol/cm^3_pw
#   • FeOOH, FeS (particulate solids): mmol/cm^3_solid
#   • pH: dimensionless free scale
species_mapping = {
    "z": ["depth_in_cm", "cm", "m"],
    "SO4": ["SO4", "mmol/cm^3_pw", "mmol/L"],
    "H2S": ["TH2S", "mmol/cm^3_pw", "mmol/L"],
    "TS2": ["TH2S", "mmol/cm^3_pw", "mmol/L"],
    "Fe2_dis": ["TFe", "mmol/cm^3_pw", "umol/L"],
    "Fe2_liq": ["Fe_aq", "mmol/cm^3_pw", "nmol/L"],
    "FeS": ["FeS", "mmol/cm^3_solid", "wt%"],
    "FeOOH": ["FeOOH", "mmol/cm^3_solid", "wt%"],
    "pH": ["pH", "pH", "pH"],
}

# ----  color map and plot style setup  --- #
colormap = "default"
colormap_override = {
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True, "zorder": -10},
}
# ----  concentration data panel  --- #
liquid_species_1 = liquid_species({
    "scale": "lin",
    "species": ["SO4"],
    # "lim": (0, 30),
    "units": "mmol/l",
    "position": "left",
})

liquid_species_2 = liquid_species({
    "scale": "lin",
    "species": ["TS2"],
    # "lim": (0, 30),
    "units": "mmol/l",
    "position": "right",
})

liquid_species_3 = liquid_species({
    "scale": "lin",
    "species": ["Fe2_dis"],
    # "lim": (0, 30),
    "units": "umol/l",
    "position": "right",
})

liquid_species_4 = liquid_species({
    "scale": "lin",
    "species": ["Fe2_liq"],
    # "lim": (0, 30),
    "units": "nmol/l",
    "position": "right",
})

solid_species = solid_species([
    {
        "scale": "lin",
        "units": "wt%",
        "relative_to": "Fe",
        "species": ["FeS"],
        "lim": (0, 2),
    },
    {
        "scale": "lin",
        "units": "wt%",
        "relative_to": "Fe",
        "species": ["FeOOH"],
        # "lim": (0, 0.12),
    },
    {
        "scale": "lin",
        "species": ["pH"],
        # "lim": (6.5, 9),
        "units": "pH",
    },
])

plot_panels = [
    [liquid_species_1, liquid_species_2, liquid_species_3,  liquid_species_4],
    [solid_species],
]
