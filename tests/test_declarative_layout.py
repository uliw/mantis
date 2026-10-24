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


