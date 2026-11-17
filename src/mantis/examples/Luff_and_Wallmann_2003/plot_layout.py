"""plot_layout.py

Declarative layout configuration for the Luff & Wallmann (2003) Hydrate Ridge model.
"""

# ---- 1. Domain & Scale Configuration ---- #
z_scale = "lin"
z_lim = (0, 0.15)  # 0 to 15 cm sediment column (in meters)

# ---- 2. Optional Comparison Against Measured Field Data ---- #
comparison_data = "data_so143_55_2.csv"
measured_data_as = "stars"  # Options: "stars", "lines", "dotted_lines", "points", "circles", etc.

# Species mapping between Mantis simulation and external measured data:
# Format: "mantis_species": [external_column_name, source_unit, target_display_unit]
species_mapping = {
    "z": ["depth_cm", "cm", "m"],
    "SO4": ["SO4_mM", "mmol/L", "mmol/L"],
    "Ca2": ["Ca_mM", "mmol/L", "mmol/L"],
    "TS2": ["H2S_mM", "mmol/L", "mmol/L"],
    "TA": ["TA_mM", "mmol/L", "mmol/L"],
    "TPO4": ["PO4_uM", "umol/L", "umol/L"],
    "TNH4": ["NH4_uM", "umol/L", "umol/L"],
    "pH": ["pH", "pH", "pH"],
    "CaCO3": ["CaCO3_wt_pct", "wt%", "wt%"],
}

# ---- 3. Colormap and Line Style Setup ---- #
colormap = "default"
colormap_override = {
    "SO4": {"color": "#1f77b4", "lw": 2},
    "CH4": {"color": "#ff7f0e", "lw": 2},
    "Ca2": {"color": "#2ca02c", "lw": 2},
    "TA": {"color": "#9467bd", "lw": 2},
    "TS2": {"color": "#8c564b", "lw": 2},
    "TPO4": {"color": "#d62728", "lw": 2},
    "TNH4": {"color": "#17becf", "lw": 2},
    "CaCO3": {"color": "#800080", "lw": 2.5},
    "Calcite": {"color": "#e377c2", "lw": 1.5, "linestyle": "--"},
    "Aragonite": {"color": "#7f7f7f", "lw": 1.5, "linestyle": "-."},
    "pH": {"color": "#bcbd22", "lw": 2},
}

# ---- 4. Panel Specifications ---- #

# Panel 1: Major porewater redox couples and calcium (SO4, CH4, Ca2)
panel_major_solutes = liquid_species({
    "scale": "lin",
    "species": ["SO4", "CH4", "Ca2"],
    "lim": (0, 75),
    "units": "mmol/l",
    "position": "left",
})

# Panel 2: Carbonate alkalinity and sulfide (TA, TS2)
panel_ta_ts2 = liquid_species({
    "scale": "lin",
    "species": ["TA", "TS2"],
    "lim": (0, 35),
    "units": "mmol/l",
    "position": "left",
})

# Panel 3: Dissolved nutrients (PO4 and NH4)
panel_nutrients_po4 = liquid_species({
    "scale": "lin",
    "species": ["TPO4"],
    "lim": (0, 60),
    "units": "umol/l",
    "position": "left",
})

panel_nutrients_nh4 = liquid_species({
    "scale": "lin",
    "species": ["TNH4"],
    "lim": (80, 240),
    "units": "umol/l",
    "position": "right",
})

# Panel 4: Carbonates (CaCO3, Calcite, Aragonite in wt%) and porewater pH
panel_carbonates = solid_species([
    {
        "scale": "lin",
        "units": ["wt"],
        "relative_to": "CaCO3",
        "species": ["CaCO3", "Calcite", "Aragonite"],
        "lim": (0, 30),
        "position": "left",
    },
    {
        "scale": "lin",
        "units": "pH",
        "species": ["pH"],
        "lim": (6.5, 9.0),
        "position": "right",
    },
])

# ---- 5. Multi-panel Layout ---- #
plot_panels = [
    [panel_major_solutes],
    [panel_ta_ts2],
    [panel_nutrients_po4, panel_nutrients_nh4],
    panel_carbonates,
]
