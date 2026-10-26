"""Unit tests for mantis_create_examples utility."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from mantis.mantis_create_examples import (
    _is_testing_environment,
    copy_examples,
    get_examples_source_path,
    main,
)


def test_is_testing_environment_detection():
    """Verify testing environment is detected when running pytest."""
    assert _is_testing_environment() is True


def test_main_skips_during_testing(capsys):
    """Verify that calling main during tests exits gracefully without copying."""
    exit_code = main([])
    assert exit_code == 0
    captured = capsys.readouterr()
    assert "execution skipped in testing environment" in captured.out


def test_get_examples_source_path_resolves():
    """Verify get_examples_source_path finds the examples directory."""
    path = get_examples_source_path()
    assert path.is_dir()
    assert (path / "velde_et_al_2016").is_dir()


def test_copy_examples_into_target_dir(tmp_path: Path):
    """Verify copying examples creates the examples tree in target directory."""
    copied_dir = copy_examples(target_dir=tmp_path)
    assert copied_dir.is_dir()
    assert (copied_dir / "velde_et_al_2016").is_dir()
    assert (copied_dir / "velde_et_al_2016" / "with_isotopes").is_dir()
    assert (
        copied_dir / "velde_et_al_2016" / "with_isotopes" / "run_velde_fast.py"
    ).is_file()

    # Verify cache and temporary files were ignored
    assert not (
        copied_dir / "velde_et_al_2016" / "with_isotopes" / "__pycache__"
    ).exists()


def test_copy_examples_overwrite_flag(tmp_path: Path):
    """Verify copy_examples handles overwrite properly."""
    copy_examples(target_dir=tmp_path)
    # Second copy without overwrite should raise FileExistsError
    with pytest.raises(FileExistsError):
        copy_examples(target_dir=tmp_path, overwrite=False)

    # Copy with overwrite=True should succeed
    dest = copy_examples(target_dir=tmp_path, overwrite=True)
    assert dest.is_dir()


def test_cli_execution_in_clean_subprocess(tmp_path: Path):
    """Verify the executable / CLI entry point executes cleanly outside test runner."""
    # Run in a clean subprocess where PYTEST_CURRENT_TEST and pytest module are absent
    clean_env = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith("PYTEST") and k != "CI"
    }

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mantis.mantis_create_examples",
            "--target-dir",
            str(tmp_path),
        ],
        cwd=str(tmp_path),
        env=clean_env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "Successfully copied MANTIS examples" in result.stdout
    assert (tmp_path / "examples" / "velde_et_al_2016").is_dir()
