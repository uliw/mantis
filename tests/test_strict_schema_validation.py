"""test_strict_schema_validation.py

Tests for the strict input file schema validation without backwards compatibility.
Verifies that any legacy formats trigger an immediate fatal exit with a helpful error message,
and that the new format is properly parsed and synthesized.
"""

import pytest
import sys
from mantis.species_parser import (
    validate_isotope_standards,
    parse_and_validate_species_dict,
    parse_boundary_condition,
    validate_and_normalize_speciation,
)


def test_legacy_isotope_standards_single_letter_fails(capsys):
    """Legacy standards keyed by 'S' without required fields must exit."""
    legacy_stds = {
        "S": {
            "standard": "VCDT",
            "ratio": 0.044162589,
            "unit": "mUr",
        }
    }
    with pytest.raises(SystemExit):
        validate_isotope_standards(legacy_stds)
    captured = capsys.readouterr()
    assert "[MANTIS ERROR]" in captured.err
    assert "Legacy isotope standard key 'S' detected" in captured.err


def test_valid_isotope_standards_passes():
    """Modern standards keyed by delta notation succeed."""
    new_stds = {
        "d34S": {
            "standard": "VCDT",
            "ratio": 0.044162589,
            "unit": "mUr",
            "element": "S",
        },
        "d18O": {
            "standard": "VSMOW",
            "ratio": 0.0020052,
            "unit": "mUr",
            "element": "O",
        },
    }
    val = validate_isotope_standards(new_stds)
    assert "d34S" in val
    assert "d18O" in val
    assert val["d34S"]["ratio"] == 0.044162589


def test_legacy_bc_key_fails(capsys):
    """Species using legacy 'bc' key must immediately exit."""
    stds = {
        "d34S": {
            "standard": "VCDT",
            "ratio": 0.044162589,
            "unit": "mUr",
            "element": "S",
        }
    }
    species = {
        "SO4": {
            "phase": "liquid",
            "include": True,
            "formula": "SO4",
            "bc": {"top": 28.2, "type": "dissolved"},
        }
    }
    with pytest.raises(SystemExit):
        parse_and_validate_species_dict(species, stds)
    captured = capsys.readouterr()
    assert "[MANTIS ERROR]" in captured.err
    assert "uses legacy 'bc' syntax" in captured.err


def test_legacy_boundary_types_dissolved_particulate_fails(capsys):
    """Using 'dissolved' or 'particulate' inside bc_top must exit."""
    stds = {
        "d34S": {
            "standard": "VCDT",
            "ratio": 0.044162589,
            "unit": "mUr",
            "element": "S",
        }
    }
    species = {
        "SO4": {
            "phase": "liquid",
            "include": True,
            "formula": "SO4",
            "bc_top": {"type": "dissolved", "top": 28.2},
        }
    }
    with pytest.raises(SystemExit):
        parse_and_validate_species_dict(species, stds)
    captured = capsys.readouterr()
    assert "[MANTIS ERROR]" in captured.err
    assert "legacy boundary type syntax" in captured.err


def test_legacy_isotope_keys_fail(capsys):
    """Using legacy 'isotope', 'delta_name', or 'isotope_of' must exit."""
    stds = {
        "d34S": {
            "standard": "VCDT",
            "ratio": 0.044162589,
            "unit": "mUr",
            "element": "S",
        }
    }
    species = {
        "SO4": {
            "phase": "liquid",
            "include": True,
            "formula": "SO4",
            "bc_top": {"Dirichlet": 28.2},
            "isotope": "SO4_32",
        }
    }
    with pytest.raises(SystemExit):
        parse_and_validate_species_dict(species, stds)
    captured = capsys.readouterr()
    assert "[MANTIS ERROR]" in captured.err
    assert "defines deprecated legacy isotope keys" in captured.err


def test_canonical_synthesis_and_bc_derivation():
    """Valid new format synthesizes companions and derives correct boundary values."""
    stds = {
        "d34S": {
            "standard": "VCDT",
            "ratio": 0.044162589,
            "unit": "mUr",
            "element": "S",
        },
        "d18O": {
            "standard": "VSMOW",
            "ratio": 0.0020052,
            "unit": "mUr",
            "element": "O",
        },
    }
    species = {
        "SO4": {
            "phase": "liquid",
            "include": True,
            "formula": "SO4",
            "diffusion": {"m0": 4.88, "m1": 0.232},
            "isotopes": {"SO4_S32": "S", "SO4_O16": "d18O"},
            "bc_top": {"Dirichlet": 28.2, "d34S": 21.0, "d18O": 0.0},
        },
        "FeS": {
            "phase": "solid",
            "include": True,
            "formula": "FeS",
            "isotopes": {"FeS_S32": "S"},
            "bc_top": {"Robin": 0.0, "d34S": 0.0},
        },
    }
    defs, bc_map, val_stds = parse_and_validate_species_dict(species, stds)

    # 1. Synthesized companions exist
    assert "SO4_S32" in defs
    assert "SO4_O16" in defs
    assert "FeS_S32" in defs

    # 2. Inherited properties
    assert defs["SO4_S32"]["phase"] == "liquid"
    assert defs["SO4_S32"]["isotope_of"] == "SO4"
    assert defs["SO4_S32"]["isotope_standard"] == "d34S"
    assert defs["SO4_S32"]["element"] == "S"
    assert defs["SO4_O16"]["isotope_standard"] == "d18O"
    assert defs["SO4_O16"]["element"] == "O"

    # 3. Companion boundary conditions derived
    assert "SO4_S32" in bc_map
    assert bc_map["SO4_S32"]["type"] == "Dirichlet"
    assert bc_map["SO4_S32"]["top"] > 0.0
    assert bc_map["SO4_S32"]["top"] < 28.2  # light mass fraction

    assert "FeS_S32" in bc_map
    assert bc_map["FeS_S32"]["type"] == "Robin"
    assert bc_map["FeS_S32"]["top"] == 0.0


def test_validate_and_normalize_speciation_explicit():
    """Verify explicit declarative speciation block validation."""
    species = {
        "TS2": {"phase": "liquid", "include": True, "bc_top": {"Dirichlet": 0.0}},
        "Fe2_total": {"phase": "mixed", "include": True, "bc_top": {"Dirichlet": 0.0}},
    }
    raw_speciation = {
        "carbonate": {
            "type": "carbonate_system",
            "primary": {"DIC": "DIC", "TA": "TA"},
            "species": ["CO2", "HCO3", "CO3", "Hplus", "pH"],
        },
        "sulfide": {
            "type": "acid_base",
            "primary": "TS2",
            "species": {
                "HS": {"fraction": "Ka1 / (Hplus + Ka1)", "include_in_TA": 1.0},
                "H2S": {"fraction": "Hplus / (Hplus + Ka1)"},
            },
            "constants": {"pKa1": 7.0},
        },
        "iron_partition": {
            "type": "partition",
            "primary": "Fe2_total",
            "fractions": {
                "dissolved": "1.0 / (1.0 + 0.001)",
                "sorbed": "0.001 / (1.0 + 0.001)",
            },
            "species": {
                "Fe2_dis": {"phase": "liquid", "fraction_of": "dissolved"},
                "Fe2_sorb": {"phase": "solid", "fraction_of": "sorbed"},
            },
        },
    }

    norm = validate_and_normalize_speciation(raw_speciation, species)
    assert "carbonate" in norm
    assert norm["carbonate"]["type"] == "carbonate_system"
    assert "sulfide" in norm
    assert norm["sulfide"]["type"] == "acid_base"
    assert "HS" in norm["sulfide"]["species"]
    assert norm["sulfide"]["species"]["HS"]["include_in_TA"] == 1.0

    # Sub-species are registered into species dictionary
    assert "HS" in species
    assert species["HS"]["include"] is False
    assert "CO3" in species
    assert species["CO3"]["include"] is False


def test_legacy_speciation_keys_rejected(capsys):
    """Verify that legacy speciation keys in species dictionary cause an immediate error."""
    from mantis.species_parser import parse_and_validate_species_dict
    species = {
        "TS2": {"phase": "liquid", "include": True, "bc_top": {"Dirichlet": 0.0}},
        "HS": {
            "phase": "liquid",
            "include": False,
            "equilibrium_of": "TS2",
            "bc_top": {"Dirichlet": 0.0},
        },
    }
    with pytest.raises(SystemExit):
        parse_and_validate_species_dict(species, {})
    captured = capsys.readouterr()
    assert "[MANTIS ERROR]" in captured.err
    assert "declares legacy speciation keys" in captured.err


def test_invalid_speciation_type_fails(capsys):
    """Verify invalid speciation type causes an immediate exit with informative error."""
    species = {"TS2": {"phase": "liquid", "include": True, "bc_top": {"Dirichlet": 0.0}}}
    invalid_spec = {
        "bad_system": {
            "type": "quantum_field",
        }
    }
    with pytest.raises(SystemExit):
        validate_and_normalize_speciation(invalid_spec, species)
    captured = capsys.readouterr()
    assert "[MANTIS ERROR]" in captured.err
    assert "invalid type 'quantum_field'" in captured.err

