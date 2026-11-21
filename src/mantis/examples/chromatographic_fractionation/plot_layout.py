"""plot_layout.py

Plot layout for the chromatographic fractionation model.
Visualizes Fe2 concentration front and chromatographic delta56Fe isotope excursion
against the analytical Ogata-Banks retarded advection-dispersion solution.
"""

# ---- 1. Domain & Scale Configuration ---- #
z_scale = "lin"
z_lim = (0, 0.20)  # 0 to 20 cm sediment column (in meters)

# ---- 2. Analytical Comparison Overlay ---- #
comparison_data = "analytical_ogata_banks.csv"
comparison_data_as = "dotted_lines"

species_mapping = {
    "z": ["z", 1.0],
    "Fe2_total": ["Fe2_total", "mmol/L", "mmol/L"],
    "Fe2_liq": ["Fe2_liq", "mmol/L", "mmol/L"],
    "Fe2_sorbed": ["Fe2_sorbed", "mmol/L", "mmol/L"],
    "d56Fe_Fe2_total": ["d56Fe_Fe2_total", 1.0],
    "d56Fe_Fe2_liq": ["d56Fe_Fe2_liq", 1.0],
    "d56Fe_Fe2_sorbed": ["d56Fe_Fe2_sorbed", 1.0],
}

# ---- 3. Color Map & Line Style Setup ---- #
colormap = "default"
colormap_override = {
    "Fe2_total": {"color": "#1f77b4", "lw": 2},
    "Fe2_liq": {"color": "#ff7f0e", "lw": 2},
    "Fe2_sorbed": {"color": "#2ca02c", "lw": 2, "linestyle": "--"},
    "d56Fe_Fe2_total": {"color": "#1f77b4", "lw": 2},
    "d56Fe_Fe2_liq": {"color": "#ff7f0e", "lw": 2},
    "d56Fe_Fe2_sorbed": {"color": "#2ca02c", "lw": 2},
}

# ---- 4. Concentration Data Panel ---- #
concentration_panel = {
    "scale": "lin",
    "species": [
        "Fe2_total",
        "Fe2_liq",
        "Fe2_sorbed",
    ],
    "units": "mmol/l",
}

# ---- 5. Isotope Data Panel ---- #
isotope_panel = {
    "species": [
        "d56Fe_Fe2_total",
        "d56Fe_Fe2_liq",
        "d56Fe_Fe2_sorbed",
    ],
    "units": "mUr",
}

plot_panels = [
    [concentration_panel],
    [isotope_panel],
]
