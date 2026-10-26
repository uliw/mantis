"""Copy bundled MANTIS examples into the current working directory."""

from __future__ import annotations

import argparse
import importlib.resources as pkg_resources
import os
import shutil
import sys
from pathlib import Path


def _is_testing_environment() -> bool:
    """Return True if running under pytest or common CI test environments."""
    if "PYTEST_CURRENT_TEST" in os.environ or "pytest" in sys.modules:
        return True
    return bool(
        os.environ.get("CI", "").lower() in ("1", "true", "yes")
        and "TEST" in "".join(os.environ.keys()).upper()
    )


def get_examples_source_path() -> Path:
    """Resolve the directory containing example files.

    First checks the installed package directory via ``importlib.resources``.
    If not found (e.g. editable install without package-data build step),
    falls back to searching repository root paths relative to the package.
    """
    import mantis

    # 1. Check if examples are packaged inside the mantis module
    try:
        pkg_root = pkg_resources.files(mantis)
        # pkg_root could be a Traversable or Path
        candidate = Path(str(pkg_root)) / "examples"
        if candidate.is_dir():
            return candidate
    except Exception:
        pass

    # 2. Check repository root relative to mantis module file
    if hasattr(mantis, "__file__") and mantis.__file__:
        mantis_dir = Path(mantis.__file__).resolve().parent
        # In a repo layout: <repo_root>/src/mantis/__init__.py
        for parent in (mantis_dir.parent.parent, mantis_dir.parent):
            candidate = parent / "examples"
            if candidate.is_dir():
                return candidate

    # 3. Fallback to cwd/examples if it exists
    cwd_candidate = Path.cwd() / "examples"
    if cwd_candidate.is_dir():
        return cwd_candidate

    raise FileNotFoundError(
        "Could not locate MANTIS examples directory in the installed package or repository."
    )


def copy_examples(target_dir: Path | str | None = None, overwrite: bool = False) -> Path:
    """Copy the examples directory into target_dir (defaults to current working directory).

    Parameters
    ----------
    target_dir : Path or str or None, optional
        Target directory where the `examples` folder should be copied into.
        Defaults to current working directory (`Path.cwd()`).
    overwrite : bool, default False
        Whether to overwrite existing files in destination directory.

    Returns
    -------
    Path
        Destination path of the copied `examples` directory.

    """
    dest_base = Path(target_dir) if target_dir else Path.cwd()
    dest_dir = dest_base / "examples"

    source_dir = get_examples_source_path()

    if dest_dir.resolve() == source_dir.resolve():
        print(
            f"Current directory already points to the source examples at: {source_dir}. "
            "No copy needed."
        )
        return dest_dir

    def _ignore_patterns(path: str, names: list[str]) -> set[str]:
        # Ignore cache files, temporary run outputs, lock files and binary state files if present
        ignored = set()
        for name in names:
            if (
                name in ("__pycache__", "TAGS", ".pytest_cache")
                or name.endswith((".pyc", ".pyo"))
                or name.startswith(".#")
                or name.startswith("#")
            ):
                ignored.add(name)
        return ignored

    shutil.copytree(
        source_dir,
        dest_dir,
        dirs_exist_ok=overwrite,
        ignore=_ignore_patterns,
    )
    return dest_dir


def main(argv: list[str] | None = None) -> int:
    """Entry point for the mantis_create_examples command line tool."""
    if _is_testing_environment():
        print("mantis_create_examples: execution skipped in testing environment.")
        return 0

    parser = argparse.ArgumentParser(
        prog="mantis_create_examples",
        description="Copy MANTIS example models and templates into current working directory.",
    )
    parser.add_argument(
        "--target-dir",
        type=Path,
        default=Path.cwd(),
        help="Destination directory to copy 'examples' into (default: current working directory).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite files in destination if 'examples' already exists.",
    )

    args = parser.parse_args(argv)

    try:
        dest_dir = copy_examples(target_dir=args.target_dir, overwrite=args.overwrite)
        print(f"Successfully copied MANTIS examples to: {dest_dir}")
        return 0
    except FileExistsError:
        print(
            f"Error: Target directory '{args.target_dir / 'examples'}' already exists.\n"
            "Use --overwrite to allow writing into the existing directory.",
            file=sys.stderr,
        )
        return 1
    except Exception as exc:
        print(f"Error copying examples: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
