# MANTIS: A Python framework to solve extremely stiff isotope enabled reaction-transport equations using a *M*onolithic *A*ssembly of *N*onlinear *T*ransport *I*nteractions *S*olver (MANTIS)

[![Build Status](https://github.com/uliw/mantis/actions/workflows/tests.yml/badge.svg)](https://github.com/uliw/mantis/actions/workflows/tests.yml)
[![PyPI Version](https://img.shields.io/pypi/v/mantis-rtm.svg)](https://pypi.org/project/mantis-rtm/)
[![Python Versions](https://img.shields.io/badge/python-3.12%20%7C%203.13%20%7C%203.14-blue.svg)](https://pypi.org/project/mantis-rtm/)
[![License: AGPL-3.0-or-later](https://img.shields.io/badge/License-AGPL_3.0--or--later-blue.svg)](https://github.com/uliw/mantis/blob/main/LICENSE)

**MANTIS** is a high-performance 1D reactive-transport modeling (RTM) framework written in Python. Unlike most reaction-transport codes, Mantis uses a fully coupled monolithic equation system. This approach allows to solve  stiff, tightly coupled, non-linear advection-diffusion-reaction networks efficiently with minimal computational overhead. Mantis supports non-uniform grids, speciation, and adsorption as well as kinetic and equilibrium isotope effects.

MANTIS models are defined through a declarative interface that is easy to read, and has an  integrated plotting system that supports asynchronous time-lapse video generation.


## How is MANTIS Different from Other RTM Packages?
Many (if not most) geochemical reactive-transport modeling packages rely on **operator splitting** (such as Sequential Non-Iterative Approach [SNIA] or Sequential Iterative Approach [SIA]) or **decoupled Picard (Gauss-Seidel) iterations**:

1. **Monolithic Coupling vs. Operator Splitting**:
   - *Operator splitting* alternates between solving pure transport and pure batch chemistry over a time step. While simple, it introduces splitting errors and often requires tiny time steps ($\Delta t$ on the order of seconds or hours) to avoid artificial oscillations or numerical divergence when reactions are fast and stiff.
   - *Picard iterations* solve each chemical species sequentially, lagging cross-species couplings from the previous iteration. At sharp reaction fronts (e.g., where sulfide meets oxygen or ferric iron), this decoupling frequently causes severe numerical lag oscillations.
   - **MANTIS solves the entire system monolithically**: Transport and reactions for all species are assembled into a single unified block-tridiagonal linear system at each Newton step. This achieves true simultaneous coupling and quadratic Newton convergence.
   
4. **Unconditional Numerical Stability**:
   - Thanks to its fully implicit monolithic formulation, MANTIS does not suffer from Courant-Friedrichs-Lewy (CFL) or diffusion time step stability restrictions, enabling long-term simulations (centuries to millennia) with large time steps.

## Current Limitations
While MANTIS provides robust and unconditionally stable solutions for stiff reactive transport, the current release has the following limitations:

1. **Time-Invariant Spatial Mesh**:
   - The spatial grid $z$ is static throughout the simulation. Dynamic grid adaptation, moving boundaries, or expanding domains (e.g., physical sediment column growth over time) are not currently supported.
2. **Time-Invariant Porosity**:
   - The porosity profile $\phi(z)$ is stationary in time (steady-state compaction). Transient sediment compaction or compaction-driven dynamic fluid expulsion over time is not currently modeled.
3. **Time-Invariant Temperature**:
   - The temperature profile of $T(z)$ is steady-state. Dynamic thermal evolution and time-dependent diffusion coefficients or rate constants driven by non-isothermal thermal transients are not currently supported.


## Installation
MANTIS requires Python 3.12 or newer. Clone the repository and install in editable mode:

```bash
pip install mantis-rtm 
mantis_create_examples.py 
```

Dependencies include `numpy`, `pandas`, `matplotlib`, `numba`, `sympy`, `chempy`, and `watchdog`.

## Quick Start Example
A quick-start tutorial is available in the `./docs` folder

## Running Tests
Run the unit test suite using `pytest`:
```bash
pytest
```

## License
This project is licensed under the GPL-3.0 License. See the [LICENSE](LICENSE) file for details.
