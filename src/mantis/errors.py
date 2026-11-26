"""Mantis-RTM: Custom exception hierarchy and clean error reporting.

Copyright (C), 2026 - Ulrich G. Wortmann

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Sequence

_BANNER = "# --------------------------------------------------------------- #"


class MantisError(Exception):
    """Base class for user-facing Mantis errors.

    By default, uncaught ``MantisError`` instances are formatted cleanly
    without a Python stack trace via :func:`mantis_excepthook`. Set
    ``MANTIS_DEBUG=1`` in the environment to force full tracebacks.
    Standard Python exceptions (not inheriting from ``MantisError``) always
    display the full stack trace.
    """

    show_traceback: bool = False

    def __init__(
        self,
        message: str,
        *,
        file_path: str | Path | None = None,
        line: int | None = None,
        context: str | None = None,
        details: Sequence[str] | None = None,
    ) -> None:
        self.message = message
        self.file_path = Path(file_path) if file_path is not None else None
        self.line = line
        self.context = context
        self.details = list(details) if details is not None else None
        super().__init__(self.format_plain())

    def _header(self, *, full_path: bool = False) -> str:
        parts: list[str] = []
        if self.context:
            parts.append(self.context)
        if self.line is not None:
            parts.append(f"line {self.line}")
        if self.file_path is not None:
            parts.append(f"in {self.file_path if full_path else self.file_path.name}")
        return " ".join(parts)

    def format_plain(self) -> str:
        """Return a concise single- or multi-line plain description (for str(exc))."""
        header = self._header(full_path=True)
        if self.details:
            body = "\n".join(self.details)
            if self.message:
                return f"{self.message}\n{body}"
            return body
        if header:
            return f"{header}: {self.message}"
        return self.message

    def format_message(self) -> str:
        """Return the banner-framed user-facing error message."""
        lines: list[str] = [_BANNER, ""]
        if self.details:
            for idx, item in enumerate(self.details):
                if idx > 0:
                    lines.append("")
                for subline in item.splitlines():
                    lines.append(f" {subline}")
        else:
            header = self._header()
            if header:
                lines.append(f" {header}:")
            for subline in self.message.splitlines():
                lines.append(f" {subline}")
        lines.extend(["", _BANNER])
        return "\n".join(lines)


class InputFileError(MantisError, ValueError):
    """Error in a user input file (invalid schema, malformed definitions, etc.)."""


class ValidationError(InputFileError):
    """Raised when chemical reaction or model specifications fail validation checks."""


class MissingDefinitionError(ValidationError):
    """A referenced species, rate constant, limiter, or function is missing from user input files."""


class ReactionParseError(MantisError, ValueError):
    """Error parsing a chemical reaction equation string or stoichiometry."""


class UnitParseError(MantisError, ValueError):
    """Error parsing or converting physical units in user input files."""


class FilePathError(MantisError, FileNotFoundError):
    """A required user input or data file could not be found at the specified path."""


class SolverError(MantisError, RuntimeError):
    """Numerical solver failed to converge or encountered an unrecoverable state."""


_original_excepthook = sys.excepthook


def mantis_excepthook(exc_type, exc_value, exc_tb) -> None:
    """Display clean banner messages for MantisError; full tracebacks for all other errors."""
    if isinstance(exc_value, MantisError) and not exc_value.show_traceback:
        if os.environ.get("MANTIS_DEBUG", "").lower() in ("1", "true", "yes"):
            _original_excepthook(exc_type, exc_value, exc_tb)
            return
        print(f"\n{exc_value.format_message()}\n", file=sys.stderr)
        return
    _original_excepthook(exc_type, exc_value, exc_tb)


def install_excepthook() -> None:
    """Install the Mantis exception hook on ``sys.excepthook``."""
    if sys.excepthook is not mantis_excepthook:
        sys.excepthook = mantis_excepthook


def uninstall_excepthook() -> None:
    """Restore the original ``sys.excepthook``."""
    sys.excepthook = _original_excepthook


__all__ = [
    "MantisError",
    "InputFileError",
    "MissingDefinitionError",
    "ValidationError",
    "ReactionParseError",
    "UnitParseError",
    "FilePathError",
    "SolverError",
    "mantis_excepthook",
    "install_excepthook",
    "uninstall_excepthook",
]
