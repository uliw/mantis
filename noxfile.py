#!/usr/bin/env -S uv run --script

# /// script
# dependencies = ["nox>=2025.10.16"]
# ///

"""Nox sessions for MANTIS."""

from __future__ import annotations

import nox

nox.needs_version = ">=2025.10.16"
nox.options.default_venv_backend = "uv|virtualenv"

PYPROJECT = nox.project.load_toml("pyproject.toml")


@nox.session(default=True)
def tests(session: nox.Session) -> None:
    """Run the unit and benchmark test suite."""
    session.install("-e", ".", *nox.project.dependency_groups(PYPROJECT, "test"))
    session.run("pytest", *session.posargs)


@nox.session(default=False)
def lint(session: nox.Session) -> None:
    """Run static type checking and linting."""
    session.install("-e", ".", *nox.project.dependency_groups(PYPROJECT, "lint"))
    session.run("mypy")
    session.run("ruff", "check", "--no-fix", *session.posargs)


@nox.session(default=False)
def docs(session: nox.Session) -> None:
    """Build the Sphinx documentation."""
    session.install("-e", ".", *nox.project.dependency_groups(PYPROJECT, "doc"))
    session.run("sphinx-build", "-b", "html", "docs", "docs/_build/html", *session.posargs)


if __name__ == "__main__":
    nox.main()
