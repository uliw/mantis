"""test_diagenetic_model.py

Unit and regression tests for DiageneticModel interface in mantis.model.
"""

from pathlib import Path

import pytest

from mantis import DiageneticModel

BASE_DIR = Path(__file__).resolve().parent.parent
EXP_DIR = BASE_DIR / "examples" / "velde_et_al_2016" / "with_isotopes"


def test_diagenetic_model_experiment_name_derivation():
    model = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
    )
    # Since caller frame is test_diagenetic_model.py, stem should be test_diagenetic_model
    assert model.experiment == "test_diagenetic_model"


def test_diagenetic_model_cli_empty_arg_query(capsys):
    model = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        experiment="test_exp",
    )
    with pytest.raises(SystemExit) as exc_info:
        model.parse_cli_args(["--solver"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "--solver" in captured.out
    assert "coupled" in captured.out or "monolithic" in captured.out


def test_diagenetic_model_setup(tmp_path):
    eq_file = tmp_path / "test_setup_exp_equations.py"
    model = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        equations=eq_file,
        layout="plot_layout.py",
        experiment="test_setup_exp",
    )
    model.configure(
        grid_setup=[["0 cm", "10 cm", "1 mm", "5 mm"]],
        t_end="1 day",
        phi=0.8,
        bioirrigation={"BI0": 0.0, "depth": "0 cm"},
    )
    mp, c, k, species_transport, z, D_mol, rxn_fn = model.setup()
    assert len(species_transport) > 0
    assert hasattr(c, "SO4")
    assert hasattr(c, "Fe3")
    assert hasattr(c, "TS2")
    assert hasattr(mp, "get_limiters")
    assert callable(mp.get_limiters)
    assert hasattr(mp, "isotope_standards")
    assert "d34S" in mp.isotope_standards
    assert mp.isotope_standards["d34S"]["unit"] == "mUr"
    assert hasattr(mp, "bc_SO4_S32")
    assert mp.bc_SO4_S32 > 0.0

    # Verify save_data produces d34S_<species>
    from mantis.diff_lib import save_data
    df, _ = save_data(mp, c, k, species_transport, z, D_mol, rxn_fn)
    assert "d34S_SO4" in df.columns
    assert "d34S_TS2" in df.columns
    assert "d34S_FeS" in df.columns


def test_diagenetic_model_missing_species_file_exits(tmp_path, capsys):
    """Verify that a missing species file produces a suitable exit message and exits."""
    model = DiageneticModel(
        base_dir=tmp_path,
        species="nonexistent_species.py",
        reactions=EXP_DIR / "fes_reactions.py",
        constants=EXP_DIR / "kinetic_constants_fast.py",
    )
    with pytest.raises(SystemExit) as exc_info:
        model.setup()
    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "[MANTIS ERROR] Could not find species definitions file: 'nonexistent_species.py'" in captured.err
    assert "Current working directory:" in captured.err
    assert "Model directory:" in captured.err


def test_diagenetic_model_missing_reactions_file_exits(tmp_path, capsys):
    """Verify that a missing reactions file produces a suitable exit message and exits."""
    model = DiageneticModel(
        base_dir=tmp_path,
        species=EXP_DIR / "species.py",
        reactions="nonexistent_reactions.py",
        constants=EXP_DIR / "kinetic_constants_fast.py",
    )
    with pytest.raises(SystemExit) as exc_info:
        model.setup()
    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "[MANTIS ERROR] Could not find declarative reactions file: 'nonexistent_reactions.py'" in captured.err
    assert "Current working directory:" in captured.err


def test_diagenetic_model_missing_constants_file_exits(tmp_path, capsys):
    """Verify that a missing constants file produces a suitable exit message and exits."""
    model = DiageneticModel(
        base_dir=tmp_path,
        species=EXP_DIR / "species.py",
        reactions=EXP_DIR / "fes_reactions.py",
        constants="nonexistent_constants.py",
    )
    with pytest.raises(SystemExit) as exc_info:
        model.setup()
    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "[MANTIS ERROR] Could not find reaction constants file: 'nonexistent_constants.py'" in captured.err


def test_diagenetic_model_finds_files_in_cwd(tmp_path, monkeypatch):
    """Verify that DiageneticModel finds simulation files in the current working directory."""
    import shutil
    shutil.copy(EXP_DIR / "species.py", tmp_path / "species.py")
    shutil.copy(EXP_DIR / "fes_reactions.py", tmp_path / "fes_reactions.py")
    shutil.copy(EXP_DIR / "kinetic_constants_fast.py", tmp_path / "kinetic_constants_fast.py")

    monkeypatch.chdir(tmp_path)
    model = DiageneticModel(
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
    )
    assert model.species_file == (tmp_path / "species.py").resolve()
    assert model.reactions_file == (tmp_path / "fes_reactions.py").resolve()
    assert model.constants_file == (tmp_path / "kinetic_constants_fast.py").resolve()


def test_diagenetic_model_finds_files_in_model_directory(tmp_path):
    """Verify that DiageneticModel finds simulation files in the directory where the model resides."""
    import shutil
    model_dir = tmp_path / "sim_subfolder"
    model_dir.mkdir()
    shutil.copy(EXP_DIR / "species.py", model_dir / "species.py")
    shutil.copy(EXP_DIR / "fes_reactions.py", model_dir / "fes_reactions.py")
    shutil.copy(EXP_DIR / "kinetic_constants_fast.py", model_dir / "kinetic_constants_fast.py")

    model = DiageneticModel(
        base_dir=model_dir,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
    )
    assert model.species_file == (model_dir / "species.py").resolve()
    assert model.reactions_file == (model_dir / "fes_reactions.py").resolve()
    assert model.constants_file == (model_dir / "kinetic_constants_fast.py").resolve()


def test_plot_summary_missing_layout_message(tmp_path, capsys):
    """Verify that plot_summary displays a suitable message when layout file is not found."""
    from mantis.diff_lib import data_container
    from mantis.model import DiageneticModelResult
    mp = data_container()
    mp.plot_name = "test_run"
    mp.max_depth = 1.0
    mp._caller_dir = str(tmp_path)
    result = DiageneticModelResult(
        mp=mp,
        c=data_container(),
        k=data_container(),
        species_list=[],
        z=[],
        D_mol=data_container(),
        diagenetic_reactions=None,
        converged="Yes",
        step=1,
        total_time=100.0,
        df=None,
    )
    res = result.plot_summary(layout_file="nonexistent_layout.py")
    assert res is None
    captured = capsys.readouterr()
    assert "[MANTIS ERROR] Plot layout file 'nonexistent_layout.py' not found." in captured.err


def test_diagenetic_model_generic_configure():
    """Verify that model.configure generically converts any parameter with physical units via Pint."""
    model = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        experiment="test_generic_cfg",
    )
    model.configure(
        custom_depth="25 cm",
        custom_duration="5 kyr",
        custom_rate="0.1 /year",
        custom_velocity="1.5 cm/yr",
        custom_diffusion="2e-5 cm^2/s",
        custom_conc="0.05 mmol/L",
        custom_flux="1.2 mol/(m^2 * year)",
        custom_tuple=("2 cm", "10 cm"),
        custom_dict={"inner_len": "5 mm", "inner_time": "2 hours"},
        plain_string="arbitrary_text",
        plain_int=42,
        plain_float=3.14,
    )
    p = model.p_dict
    assert p["custom_depth"] == pytest.approx(0.25)
    assert p["custom_duration"] == pytest.approx(5000.0 * 365.25 * 86400.0, rel=1e-3)
    assert p["custom_rate"] == pytest.approx(0.1 / (365.25 * 86400.0), rel=1e-3)
    assert p["custom_velocity"] == pytest.approx(0.015 / (365.25 * 86400.0), rel=1e-3)
    assert p["custom_diffusion"] == pytest.approx(2e-9)
    assert p["custom_conc"] == pytest.approx(0.05)
    assert p["custom_flux"] == pytest.approx(1.2 / (365.25 * 86400.0), rel=1e-3)
    assert p["custom_tuple"] == (pytest.approx(0.02), pytest.approx(0.1))
    assert p["custom_dict"]["inner_len"] == pytest.approx(0.005)
    assert p["custom_dict"]["inner_time"] == pytest.approx(7200.0)
    assert p["plain_string"] == "arbitrary_text"
    assert p["plain_int"] == 42
    assert p["plain_float"] == 3.14


def test_resolve_and_prepare_equations_triggers_regen(tmp_path):
    """Verify resolve_and_prepare_equations regenerates equations when species, reactions, constants, or force_regen is updated."""
    import time

    from mantis.generate_equation_system import resolve_and_prepare_equations

    sp_file = tmp_path / "species.py"
    rxn_file = tmp_path / "test_reactions.py"
    kc_file = tmp_path / "kinetic_constants_fast.py"
    eq_file = tmp_path / "test_equations.py"

    import shutil
    shutil.copy(EXP_DIR / "species.py", sp_file)
    shutil.copy(EXP_DIR / "fes_reactions.py", rxn_file)
    shutil.copy(EXP_DIR / "kinetic_constants_fast.py", kc_file)

    p_dict = {
        "experiment": "test",
        "species": str(sp_file),
        "reactions": str(rxn_file),
        "kinetic_constants": str(kc_file),
        "equations": str(eq_file),
    }

    # First call: generates equation file
    eq_mod, _, _ = resolve_and_prepare_equations(p_dict, base_dir=tmp_path)
    assert eq_file.exists()
    mtime_init = eq_file.stat().st_mtime

    # Immediate second call without changes: no regeneration
    time.sleep(0.05)
    eq_mod2, _, _ = resolve_and_prepare_equations(p_dict, base_dir=tmp_path)
    assert eq_file.stat().st_mtime == mtime_init

    # Touch species.py to make it newer
    time.sleep(0.05)
    sp_file.touch()
    eq_mod3, _, _ = resolve_and_prepare_equations(p_dict, base_dir=tmp_path)
    assert eq_file.stat().st_mtime > mtime_init
    mtime_after_sp = eq_file.stat().st_mtime

    # Touch kinetic_constants to make it newer
    time.sleep(0.05)
    kc_file.touch()
    eq_mod4, _, _ = resolve_and_prepare_equations(p_dict, base_dir=tmp_path)
    assert eq_file.stat().st_mtime > mtime_after_sp
    mtime_after_kc = eq_file.stat().st_mtime

    # Force regen flag
    time.sleep(0.05)
    p_dict["force_regen"] = True
    eq_mod5, _, _ = resolve_and_prepare_equations(p_dict, base_dir=tmp_path)
    assert eq_file.stat().st_mtime > mtime_after_kc


def test_exponential_bioturbation_profile():
    """Verify exponential bioturbation profile formula and configuration."""
    import numpy as np
    from mantis.diff_lib import compute_exponential_db

    z = np.array([0.0, 0.03, 0.06])
    Db0 = 10.0
    x_atten = 0.03
    db = compute_exponential_db(z, Db0, x_atten)
    np.testing.assert_allclose(db[0], 10.0)
    np.testing.assert_allclose(db[1], 10.0 / np.e)
    np.testing.assert_allclose(db[2], 10.0 / (np.e**2))


def test_all_species_exported_and_saved(tmp_path):
    """Verify all species defined in species.py are exported to CSV/df and state .npz."""
    import numpy as np
    from mantis.diff_lib import save_data, save_state

    model = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        layout="plot_layout.py",
        experiment="test_export_exp",
    )
    model.configure(
        grid_setup=[["0 cm", "10 cm", "1 mm", "5 mm"]],
        t_end="1 day",
        phi=0.8,
        bioturbation={"BT0": "10 cm^2/year", "profile": "exponential", "attenuation": "3 cm"},
    )
    mp, c, k, species_transport, z, D_mol, rxn_fn = model.setup()

    # Fe2_total initial value
    c.Fe2_total.setValue(1.0)

    df, _ = save_data(mp, c, k, species_transport, z, D_mol, rxn_fn)
    # Fe2_liq should be present in df even though include=False in species.py
    assert "c_Fe2_liq" in df.columns
    assert "c_Fe2_total" in df.columns
    assert np.all(df["c_Fe2_liq"] > 0.0)

    # State file should also contain Fe2_liq
    state_file = tmp_path / "test_state.npz"
    save_state(c, str(state_file), mp=mp)
    assert state_file.exists()
    loaded = np.load(str(state_file))
    assert "Fe2_liq" in loaded
    assert np.all(loaded["Fe2_liq"] > 0.0)


def test_diffusive_boundary_layer_stencil():
    """Verify that a diffusive boundary layer (DBL) reduces interfacial conductance."""
    import numpy as np
    from mantis.monolithic_solver import build_native_1d_transport_stencil

    z = np.linspace(0.001, 0.1, 10)
    dx = np.full(10, 0.01)
    phi = 0.8
    D = np.full(10, 1e-9)
    w = 1e-10

    # 1. Pure Dirichlet (no DBL)
    bc_dirichlet = {"type": "dissolved", "top": 28.0}
    ab_dir, b_dir = build_native_1d_transport_stencil(
        z=z, dx=dx, phi=phi, D_cell=D, w=w, bc_props=bc_dirichlet
    )

    # 2. Robin with DBL (delta = 0.5 mm = 0.0005 m)
    bc_dbl = {"type": "dissolved", "top": 28.0, "dbl": 0.0005}
    ab_dbl, b_dbl = build_native_1d_transport_stencil(
        z=z, dx=dx, phi=phi, D_cell=D, w=w, bc_props=bc_dbl
    )

    assert ab_dir.shape == (3, 10)
    assert ab_dbl.shape == (3, 10)
    # The added DBL series resistance must lower interfacial exchange conductance
    assert b_dbl[0] < b_dir[0]
    assert ab_dbl[1, 0] < ab_dir[1, 0]


def test_diffusive_boundary_layer_configure():
    """Verify that model.configure correctly parses and converts DBL units."""
    import numpy as np
    from mantis import DiageneticModel
    model = DiageneticModel()
    model.configure(diffusive_boundary_layer="0.05 cm")
    assert np.isclose(model.p_dict["diffusive_boundary_layer"], 0.0005)
    assert np.isclose(model.p_dict["dbl"], 0.0005)


def test_coupled_speciation_auto_numerical_jacobian(tmp_path):
    """Verify that coupled non-linear speciation auto-enables numerical jacobian unless explicitly set."""
    eq_file = tmp_path / "test_auto_jac_eq.py"
    model = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        equations=eq_file,
    )
    model.configure(enable_iron_speciation=True, beta_FeHS=251.2)
    model.setup()
    assert model.p_dict["use_numerical_jacobian"] is True
    assert model.mp.use_numerical_jacobian is True

    # When explicitly overridden to False by user, should respect user setting
    model2 = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        equations=eq_file,
    )
    model2.configure(enable_iron_speciation=True, beta_FeHS=251.2, use_numerical_jacobian=False)
    model2.setup()
    assert model2.p_dict["use_numerical_jacobian"] is False
    assert getattr(model2.mp, "use_numerical_jacobian", False) is False


def test_uncoupled_fixed_ph_selects_analytical_jacobian(tmp_path):
    """Verify that uncoupled / static scalar speciation with fixed pH selects analytical Jacobian."""
    eq_file = tmp_path / "test_ana_jac_eq.py"
    model = DiageneticModel(
        base_dir=EXP_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        equations=eq_file,
    )
    # Default Velde 2016 has fixed pH, constant Fe2 partition -> uncoupled
    model.setup()
    assert model.p_dict["use_numerical_jacobian"] is False
    assert model.mp.use_numerical_jacobian is False






