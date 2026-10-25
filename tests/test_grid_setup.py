"""Unit tests for multi-section mesh generation (make_grid_from_sections) and DiageneticModel integration."""

import numpy as np
import pytest
from mantis.diff_lib import make_grid_from_sections
from mantis import DiageneticModel


def test_single_uniform_section():
    """Verify single uniform section generates exact cell counts and dimensions."""
    mesh, z = make_grid_from_sections([["0 cm", "10 cm", "1 mm", "1 mm"]], verbose=False)
    assert mesh.numberOfCells == 100
    assert np.allclose(mesh.dx, 0.001)
    assert np.isclose(mesh.faceCoordinates[-1], 0.10)
    assert np.isclose(mesh.faceCoordinates[0], 0.0)
    assert len(z) == 100


def test_single_geometric_section():
    """Verify single non-uniform section expands monotonically and preserves domain length."""
    mesh, z = make_grid_from_sections([["0 cm", "50 cm", "0.1 mm", "1 cm"]], verbose=False)
    assert mesh.numberOfCells > 10
    assert np.all(mesh.dx > 0)
    assert np.isclose(mesh.faceCoordinates[0], 0.0)
    assert np.isclose(mesh.faceCoordinates[-1], 0.50)
    # Spacing should increase monotonically
    assert mesh.dx[-1] > mesh.dx[0]


def test_multisection_user_example():
    """Verify user 4-section grid setup matches boundary faces and stays stable."""
    grid_setup = [
        ["0 cm", "1 cm", "50 um", "1 mm"],
        ["1 cm", "5 cm", "1 mm", "5 mm"],
        ["5 cm", "50 cm", "5 mm", "1 cm"],
        ["50 cm", "52 cm", "1 cm", "0.1 mm"],
    ]
    mesh, z = make_grid_from_sections(grid_setup, verbose=False)
    faces = mesh.faceCoordinates

    assert np.isclose(faces[0], 0.0)
    assert np.isclose(faces[-1], 0.52)
    assert np.all(mesh.dx > 0)

    # Check exact boundary face snapping
    for boundary in [0.0, 0.01, 0.05, 0.50, 0.52]:
        assert np.any(np.isclose(faces, boundary, atol=1e-9)), f"Missing boundary {boundary} in faces"


def test_multisection_spline_mode():
    """Verify monotonic PCHIP spline interpolation mode across sections."""
    grid_setup = [
        ["0 cm", "1 cm", "50 um", "1 mm"],
        ["1 cm", "5 cm", "1 mm", "5 mm"],
        ["5 cm", "50 cm", "5 mm", "1 cm"],
    ]
    mesh, z = make_grid_from_sections(grid_setup, smoothing="spline", verbose=False)
    assert np.all(mesh.dx > 0)
    assert np.isclose(mesh.faceCoordinates[-1], 0.50)
    for boundary in [0.0, 0.01, 0.05, 0.50]:
        assert np.any(np.isclose(mesh.faceCoordinates, boundary, atol=1e-9))


def test_discontinuous_boundary_smoothing():
    """Verify that discontinuous spacing step changes are smoothly transitioned."""
    grid_setup = [
        ["0 cm", "5 cm", "1 mm", "1 mm"],
        ["5 cm", "50 cm", "5 mm", "5 mm"],
    ]
    mesh, z = make_grid_from_sections(grid_setup, smooth_discontinuities=True, verbose=False)
    assert np.isclose(mesh.faceCoordinates[0], 0.0)
    assert np.isclose(mesh.faceCoordinates[-1], 0.50)
    # Check that adjacent ratio does not jump 5x instantaneously
    ratios = mesh.dx[1:] / mesh.dx[:-1]
    assert np.max(ratios) < 3.0


def test_invalid_sections_raise():
    """Verify validation checks on invalid inputs."""
    # Empty
    with pytest.raises(ValueError, match="non-empty"):
        make_grid_from_sections([])

    # Wrong number of elements
    with pytest.raises(ValueError, match="exactly 4 elements"):
        make_grid_from_sections([["0 cm", "10 cm", "1 mm"]])

    # Stop <= start
    with pytest.raises(ValueError, match="strictly greater"):
        make_grid_from_sections([["10 cm", "5 cm", "1 mm", "1 mm"]])

    # Non-positive dz
    with pytest.raises(ValueError, match="positive"):
        make_grid_from_sections([["0 cm", "10 cm", "-1 mm", "1 mm"]])

    # Non-contiguous
    with pytest.raises(ValueError, match="Sections must be contiguous"):
        make_grid_from_sections([
            ["0 cm", "5 cm", "1 mm", "1 mm"],
            ["6 cm", "10 cm", "1 mm", "1 mm"],
        ])


def test_diagenetic_model_grid_setup(tmp_path):
    """Verify DiageneticModel automatically derives cumulative max_depth and default display_length."""
    model = DiageneticModel()
    model.configure(
        grid_setup=[
            ["0 cm", "2 cm", "0.1 mm", "0.1 mm"],
            ["2 cm", "15 cm", "0.1 mm", "2 mm"],
        ],
        phi=0.8,
    )
    # Before setup, configure should have parsed grid_setup with unit conversion
    assert len(model.mp.grid_setup) == 2
    assert np.isclose(model.mp.grid_setup[0][1], 0.02)
    assert np.isclose(model.mp.grid_setup[1][1], 0.15)


def test_diagenetic_model_display_length_default():
    """Verify display_length defaults to grid length (cumulative max_depth)."""
    grid = [
        ["0 cm", "1 cm", "0.1 mm", "0.1 mm"],
        ["1 cm", "25 cm", "0.1 mm", "1 mm"],
    ]
    mesh, z = make_grid_from_sections(grid, verbose=False)
    assert np.isclose(mesh.faceCoordinates[-1], 0.25)
