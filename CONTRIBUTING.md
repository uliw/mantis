# Contributing to MANTIS

Thank you for your interest in contributing to **MANTIS** (Monolithic Assembly of Nonlinear Transport Interactions Solver)! We welcome bug reports, documentation improvements, feature requests, and code contributions.

---

## Code of Conduct

By participating in this project, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

---

## Getting Started

### 1. Prerequisites
- **Python 3.13+**
- Recommended environment and package runner: [`uv`](https://docs.astral.sh/uv/) (or standard `pip` / `virtualenv`)
- Recommended test runner: [`nox`](https://nox.thea.codes/)

### 2. Fork and Clone
```bash
git clone https://github.com/uliw/mantis.git
cd mantis
```

### 3. Development Setup
Create and activate an environment, then install MANTIS in editable mode with development dependencies:

```bash
# Using uv (recommended)
uv venv
source .venv/bin/activate
uv pip install -e ".[dev,test,lint,doc]"

# Install pre-commit hooks
pre-commit install
```

---

## Development Workflow

We use standard toolchains configured via `nox` and `pyproject.toml`:

### Running Tests
Run the test suite using `pytest` or `nox`:
```bash
# Direct pytest run
pytest

# Or via nox
nox -s tests
```

### Code Formatting & Linting
We enforce formatting and static checks via **ruff**, **mypy**, and **pre-commit**:
```bash
# Run linting session (mypy and ruff check)
nox -s lint

# Or directly:
ruff check .
ruff format .
mypy
```

### Building Documentation
To test Sphinx documentation builds locally:
```bash
nox -s docs
# Documentation is built into docs/_build/html
```

---

## Submitting Changes

1. **Create a branch**: Branch off `staging` (or `main`) with a descriptive name:
   ```bash
   git checkout -b feature/my-new-feature
   ```
2. **Make focused commits**: Keep commits logical, focused, and accompanied by clear commit messages.
3. **Add tests**: Ensure new features or bug fixes include corresponding unit or regression tests in `tests/`.
4. **Verify checks**: Ensure tests, linters, and type checkers pass before submitting:
   ```bash
   nox -s tests lint
   ```
5. **Open a Pull Request**: Submit your pull request against the `staging` branch. Provide a concise summary of the motivation, changes made, and any related issue numbers.

---

## Reporting Issues

If you find a bug or have a feature suggestion, please open an issue on the [GitHub Issue Tracker](https://github.com/uliw/mantis/issues):
- Provide a clear, descriptive title.
- For bugs, include a minimal reproducible example (code snippet and dataset if applicable), environment details (Python version, OS, key package versions), and the complete traceback.
