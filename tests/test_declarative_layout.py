import os
import tempfile
import numpy as np
import pandas as pd
import pytest

from mantis.plot_data_new import load_layout_from_file, plot


@pytest.fixture
def sample_df():
    z = np.linspace(0, 1, 50)
    phi = 0.8 * np.ones_like(z)
    return pd.DataFrame({
        "z": z,
        "phi": phi,
        "c_SO4": np.linspace(28.0, 1.0, 50),
        "c_HS": np.linspace(0.01, 5.0, 50),
        "c_O2": np.linspace(0.2, 0.0, 50),
        "c_Fe3": np.linspace(10.0, 1.0, 50),
        "c_FeS": np.linspace(0.1, 8.0, 50),
        "c_S0": np.linspace(0.05, 0.5, 50),
        "f_SO4": -np.linspace(1e-6, 1e-8, 50),
        "f_FeS": np.linspace(1e-7, 1e-9, 50),
        "f_POC_slow": -np.linspace(1e-8, 1e-10, 50),
        "f_POC_fast": -np.linspace(1e-7, 1e-9, 50),
        "D_bio": 1e-9 * np.ones_like(z),
        "d34S_SO4": np.linspace(21.0, 45.0, 50),
        "d34S_FeS": np.linspace(-15.0, -10.0, 50),
        "d18O_SO4": np.linspace(10.0, 14.0, 50),
    })


def test_declarative_layout_basic_liquid_and_solid(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
liquid_species = {
    "scale": "log",
    "species": ["SO4", "HS"],
    "lim": (1e-4, 50),
    "units": "mmol/l",
}

solid_species = [
    {
        "scale": "lin",
        "units": ["wt"],
        "relative_to": "Fe",
        "species": ["Fe3", "FeS"],
        "lim": (0, 2),
    }
]
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        assert "first_subplot" in plt_desc
        subplot = plt_desc["first_subplot"]
        assert len(subplot["left"]) == 2
        assert len(subplot["right1"]) == 2
        assert subplot["left_ylabel"] == "Concentration [mmol/l]"
        assert subplot["right1_ylabel"] == "[wt Fe]"
        assert subplot["yscale"] == "log"
        assert subplot["ylim"] == (1e-4, 50)
        assert subplot["right1_ylim"] == (0, 2)


def test_declarative_layout_multiple_solid_axes(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
liquid_species = {
    "species": ["SO4"],
}

solid_species = [
    {
        "units": "wt%",
        "relative_to": "Fe",
        "species": ["FeS"],
    },
    {
        "units": "wt%",
        "relative_to": "S",
        "species": ["S0"],
    },
]
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        subplot = plt_desc["first_subplot"]
        assert "right1" in subplot
        assert "right2" in subplot
        assert subplot["right1_ylabel"] == "[wt% Fe]"
        assert subplot["right2_ylabel"] == "[wt% S]"


def test_declarative_layout_reaction_rates_and_poc(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
reaction_rates = {
    "species": ["SO4", "FeS", "POC", "D_bio"],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        assert "second_subplot" in plt_desc
        subplot = plt_desc["second_subplot"]
        assert len(subplot["left"]) == 4
        assert "symlog" in subplot.get("options-left", "")


def test_declarative_layout_isotopes_dual_axes_and_broken_axis(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
isotopes = {
    "species": ["d34S_SO4", "d34S_FeS", "d18O_SO4"],
    "lim": (-45, -37, 20, 30),
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        assert "third_subplot" in plt_desc
        subplot = plt_desc["third_subplot"]
        # Left has d34S (2 curves), Right1 has d18O (1 curve)
        assert len(subplot["left"]) == 2
        assert "right1" in subplot
        assert len(subplot["right1"]) == 1
        assert r"\delta^{34}" in subplot["left_ylabel"]
        assert r"\delta^{18}" in subplot["right1_ylabel"]
        assert "broken_axis(-45, -37, 20, 30)" in str(subplot.get("options-left"))


def test_declarative_layout_comparison_data_mapping(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        comp_df = pd.DataFrame({
            "depth_cm": np.linspace(0, 100, 20),
            "C_SO4": np.linspace(28.0, 2.0, 20),
            "C_FeS": np.linspace(0.01, 0.5, 20),
        })
        comp_file = os.path.join(tmpdir, "external_data.csv")
        comp_df.to_csv(comp_file, index=False)

        layout_content = f"""
comparison_data = "{comp_file}"

species_mapping = {{
    "z": ["depth_cm", 1 / 100],
    "SO4": ["C_SO4", 1.0],
    "FeS": ["C_FeS", "wt"],
}}

liquid_species = {{
    "species": ["SO4"],
}}

solid_species = [
    {{
        "species": ["FeS"],
        "relative_to": "Fe",
    }}
]
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        subplot = plt_desc["first_subplot"]
        # In left axis: SO4 (mantis) + SO4 (data) -> 2 items
        assert len(subplot["left"]) == 2
        # Comparison series with custom x has format: [ext_z, ext_y, label, style]
        assert subplot["left"][1][3]["linestyle"] == ":"
        # In right1 axis: FeS (mantis) + FeS (data) -> 2 items
        assert len(subplot["right1"]) == 2
        assert subplot["right1"][1][3]["linestyle"] == ":"


def test_declarative_layout_missing_mapped_key_raises(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        comp_df = pd.DataFrame({
            "depth_cm": [0, 10],
            "C_NONEXISTENT": [1, 2],
        })
        comp_file = os.path.join(tmpdir, "external_data.csv")
        comp_df.to_csv(comp_file, index=False)

        # Mapped external column doesn't exist
        layout_content = f"""
comparison_data = "{comp_file}"
species_mapping = {{
    "z": ["depth_cm", 0.01],
    "SO4": ["UNKNOWN_COL", 1.0],
}}
liquid_species = {{"species": ["SO4"]}}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        with pytest.raises(KeyError, match="Mapped external column 'UNKNOWN_COL' does not exist"):
            load_layout_from_file(sample_df, layout_path)


def test_declarative_layout_missing_mantis_species_raises(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        comp_df = pd.DataFrame({"depth_cm": [0, 10], "C_FOO": [1, 2]})
        comp_file = os.path.join(tmpdir, "external_data.csv")
        comp_df.to_csv(comp_file, index=False)

        # Mapped mantis species doesn't exist in simulation dataframe
        layout_content = f"""
comparison_data = "{comp_file}"
species_mapping = {{
    "z": ["depth_cm", 0.01],
    "NONEXISTENT_SPECIES": ["C_FOO", 1.0],
}}
liquid_species = {{"species": ["SO4"]}}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        with pytest.raises(KeyError, match="Species mapping key 'NONEXISTENT_SPECIES' does not exist"):
            load_layout_from_file(sample_df, layout_path)


def test_declarative_layout_colormap_override(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
colormap_override = {
    "SO4": {"color": "magenta", "lw": 4, "fill_only": True},
}

liquid_species = {
    "species": ["SO4"],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        style = plt_desc["first_subplot"]["left"][0][2]
        assert style["color"] == "magenta"
        assert style["lw"] == 4
        assert style["fill_only"] is True


def test_empty_declarative_layout_raises(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_path = os.path.join(tmpdir, "empty_layout.py")
        with open(layout_path, "w") as f:
            f.write("# empty layout\n")

        with pytest.raises(ValueError, match="must define at least one of"):
            load_layout_from_file(sample_df, layout_path)


def test_reaction_rates_fill_only_preserved(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
colormap_override = {
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True, "zorder": -10},
    "SO4": {"color": "cornflowerblue", "fillonly": True},
}
reaction_rates = {
    "species": ["SO4", "D_bio", ("FeS", {"fill_only": True, "color": "orange"})],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        subplot = plt_desc["second_subplot"]
        # SO4:
        so4_style = subplot["left"][0][2]
        assert so4_style["fill_only"] is True
        assert so4_style["color"] == "cornflowerblue"
        # D_bio:
        dbio_style = subplot["left"][1][2]
        assert dbio_style["fill_only"] is True
        assert dbio_style["color"] == "#bbbbbb"
        assert dbio_style["zorder"] == -10
        # FeS inline:
        fes_style = subplot["left"][2][2]
        assert fes_style["fill_only"] is True
        assert fes_style["color"] == "orange"


def test_isotope_resolution_flexible_naming(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Layout uses d34S_SO4, dataframe only has d_SO4
        df_d = sample_df.copy()
        df_d["d_SO4"] = df_d["d34S_SO4"]
        df_d.drop(columns=["d34S_SO4"], inplace=True)

        layout_content = """
isotopes = {
    "species": ["d34S_SO4"],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(df_d, layout_path)
        assert len(plt_desc["third_subplot"]["left"]) == 1
        # Check that it resolved to df_d["d_SO4"] data
        np.testing.assert_allclose(plt_desc["third_subplot"]["left"][0][0], df_d["d_SO4"])


def test_plot_execution_with_fill_only(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
colormap_override = {
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True, "zorder": -10},
}
liquid_species = {
    "scale": "log",
    "species": ["SO4"],
}
reaction_rates = {
    "scale": "symlog",
    "species": ["SO4", "D_bio"],
}
isotopes = {
    "species": ["d34S_SO4"],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        pdf_path = os.path.join(tmpdir, "test_output.pdf")
        plt_desc = load_layout_from_file(sample_df, layout_path)
        plot(sample_df, display_length=1.0, outfile=pdf_path, show=False, plot_description=plt_desc)
        assert os.path.exists(pdf_path)
        assert os.path.getsize(pdf_path) > 0


def test_shared_color_for_related_species_across_panels(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
liquid_species = {
    "species": ["SO4"],
}
solid_species = [
    {
        "species": ["FeS"],
        "relative_to": "Fe",
    }
]
reaction_rates = {
    "species": ["FeS", "SO4"],
}
isotopes = {
    "species": ["d34S_SO4", "d34S_FeS"],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        
        # Panel 1: SO4 (left) and FeS (right1)
        so4_p1_color = plt_desc["first_subplot"]["left"][0][2]["color"]
        fes_p1_color = plt_desc["first_subplot"]["right1"][0][2]["color"]

        # Panel 2: FeS (index 0) and SO4 (index 1)
        fes_p2_color = plt_desc["second_subplot"]["left"][0][2]["color"]
        so4_p2_color = plt_desc["second_subplot"]["left"][1][2]["color"]

        # Panel 3: d34S_SO4 (index 0) and d34S_FeS (index 1)
        d_so4_p3_color = plt_desc["third_subplot"]["left"][0][2]["color"]
        d_fes_p3_color = plt_desc["third_subplot"]["left"][1][2]["color"]

        # Assert shared color across all panels
        assert so4_p1_color == so4_p2_color == d_so4_p3_color, (
            f"SO4 colors mismatch: P1={so4_p1_color}, P2={so4_p2_color}, P3={d_so4_p3_color}"
        )
        assert fes_p1_color == fes_p2_color == d_fes_p3_color, (
            f"FeS colors mismatch: P1={fes_p1_color}, P2={fes_p2_color}, P3={d_fes_p3_color}"
        )


def test_shared_color_with_custom_colormap_cycling(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a dataframe with novel species names not in geochem defaults
        df_custom = sample_df.copy()
        df_custom["c_NovelA"] = df_custom["c_SO4"]
        df_custom["c_NovelB"] = df_custom["c_HS"]
        df_custom["f_NovelA"] = df_custom["f_SO4"]
        df_custom["f_NovelB"] = df_custom["f_FeS"]
        df_custom["d34S_NovelA"] = df_custom["d34S_SO4"]

        layout_content = """
colormap = "Set1"
liquid_species = {
    "species": ["NovelA", "NovelB"],
}
reaction_rates = {
    "species": ["NovelB", "NovelA"],
}
isotopes = {
    "species": ["d34S_NovelA"],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(df_custom, layout_path)

        a_p1_color = plt_desc["first_subplot"]["left"][0][2]["color"]
        b_p1_color = plt_desc["first_subplot"]["left"][1][2]["color"]

        b_p2_color = plt_desc["second_subplot"]["left"][0][2]["color"]
        a_p2_color = plt_desc["second_subplot"]["left"][1][2]["color"]

        d_a_p3_color = plt_desc["third_subplot"]["left"][0][2]["color"]

        assert a_p1_color == a_p2_color == d_a_p3_color
        assert b_p1_color == b_p2_color
        assert a_p1_color != b_p1_color


def test_solid_species_relative_to_none_or_missing(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        sample_df_with_ph = sample_df.copy()
        sample_df_with_ph["pH"] = 7.5

        layout_content = """
solid_species = [
    {
        "scale": "lin",
        "units": "wt",
        "relative_to": "Fe",
        "species": ["FeS"],
    },
    {
        "scale": "lin",
        "units": "pH",
        "relative_to": None,
        "species": ["pH"],
    },
    {
        "scale": "lin",
        "units": "dimensionless",
        "species": ["phi"],
    },
]
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df_with_ph, layout_path)
        first_sp = plt_desc["first_subplot"]

        # right1 has relative_to="Fe"
        assert first_sp["right1_ylabel"] == "[wt Fe]"
        assert first_sp["right1"][0][1] == "FeS [wt% Fe]"

        # right2 has relative_to=None
        assert first_sp["right2_ylabel"] == "[pH]"
        assert first_sp["right2"][0][1] == "pH"

        # right3 has missing relative_to
        assert first_sp["right3_ylabel"] == "[dimensionless]"
        assert first_sp["right3"][0][1] == "phi"


def test_custom_plot_panels_arbitrary_ordering_and_assignment(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
liquid_species_1 = liquid_species({
    "scale": "log",
    "species": ["SO4", "HS"],
    "lim": (1e-4, 30),
    "units": "mmol/l",
    "position": "left",
})

solid_species_1 = solid_species({
    "scale": "lin",
    "units": "wt%",
    "relative_to": "Fe",
    "species": ["FeS"],
    "lim": (0, 2),
    "position": "right",
})

reaction_rates_1 = reaction_rates({
    "species": ["SO4", "FeS"],
    "position": "left",
})

isotopes_1 = isotopes({
    "species": ["d34S_SO4", "d34S_FeS"],
    "lim": (-45, -37, 20, 30),
    "position": "left",
})

liquid_species_2 = liquid_species({
    "scale": "lin",
    "species": ["O2"],
    "units": "mmol/l",
})

# Flexible panel composition:
# Panel 1: Isotopes on top
# Panel 2: Liquids (left) + Solids (right)
# Panel 3: Reaction rates
# Panel 4: Another liquid panel (e.g. O2)
plot_panels = [
    [isotopes_1],
    [liquid_species_1, solid_species_1],
    [reaction_rates_1],
    [liquid_species_2],
]
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)

        assert "first_subplot" in plt_desc
        assert "second_subplot" in plt_desc
        assert "third_subplot" in plt_desc
        assert "fourth_subplot" in plt_desc

        # Panel 1 is isotopes
        p1 = plt_desc["first_subplot"]
        assert len(p1["left"]) == 2
        assert r"\delta^{34}" in p1["left_ylabel"]
        assert "broken_axis(-45, -37, 20, 30)" in str(p1.get("options-left"))

        # Panel 2 is liquid + solid
        p2 = plt_desc["second_subplot"]
        assert len(p2["left"]) == 2
        assert p2["left_ylabel"] == "Concentration [mmol/l]"
        assert "right1" in p2
        assert len(p2["right1"]) == 1
        assert p2["right1_ylabel"] == "[wt% Fe]"

        # Panel 3 is reaction rates
        p3 = plt_desc["third_subplot"]
        assert len(p3["left"]) == 2
        assert "reaction rate" in p3["left_ylabel"]

        # Panel 4 is liquid_species_2 (O2)
        p4 = plt_desc["fourth_subplot"]
        assert len(p4["left"]) == 1
        assert p4["left"][0][1] == r"O$_{2}$"
        assert p4["yscale"] == "linear"

        # Verify unified entity color sharing across the panels:
        # SO4 in p1 (isotopes), p2 (liquid), p3 (rates) must share the same color
        so4_iso_color = p1["left"][0][2]["color"]
        so4_liq_color = p2["left"][0][2]["color"]
        so4_rr_color = p3["left"][0][2]["color"]
        assert so4_iso_color == so4_liq_color == so4_rr_color

        # FeS in p1 (isotopes), p2 (solids), p3 (rates) must share the same color
        fes_iso_color = p1["left"][1][2]["color"]
        fes_sol_color = p2["right1"][0][2]["color"]
        fes_rr_color = p3["left"][1][2]["color"]
        assert fes_iso_color == fes_sol_color == fes_rr_color

        # Render PDF to confirm end-to-end plotting works without errors
        pdf_path = os.path.join(tmpdir, "panels_test.pdf")
        plot(sample_df, display_length=1.0, outfile=pdf_path, show=False, plot_description=plt_desc)
        assert os.path.exists(pdf_path)
        assert os.path.getsize(pdf_path) > 0


def test_z_scale_and_z_lim_options(sample_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        layout_content = """
z_scale = "linear"
z_lim = (0.05, 0.8)

liquid_species = {
    "species": ["SO4"],
}
"""
        layout_path = os.path.join(tmpdir, "plot_layout.py")
        with open(layout_path, "w") as f:
            f.write(layout_content)

        plt_desc = load_layout_from_file(sample_df, layout_path)
        first_sp = plt_desc["first_subplot"]
        assert first_sp["xscale"] == "linear"
        assert first_sp["xlim"] == (0.05, 0.8)

        # Verify plot application sets axes limits
        fig, axes = plot(sample_df, display_length=1.0, outfile=None, show=False, plot_description=plt_desc)
        ax = axes[0]
        assert ax.get_xscale() == "linear"
        xlim = ax.get_xlim()
        assert pytest.approx(xlim[0], rel=1e-3) == 0.05
        assert pytest.approx(xlim[1], rel=1e-3) == 0.8





