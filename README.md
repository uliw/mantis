# <img src="mantis.svg" alt="MANTIS Logo" width="20%" style="vertical-align: middle; max-width: 200px; display: inline-block; margin-right: 10px;" align="left" />  A Python framework to solve stiff reaction-transport systems in seconds using a *M*onolithic *A*ssembly of *N*onlinear *T*ransport *I*nteractions *S*olver (MANTIS)

[![Build Status](https://github.com/uliw/mantis/actions/workflows/tests.yml/badge.svg)](https://github.com/uliw/mantis/actions/workflows/tests.yml)
[![Type Checking](https://github.com/uliw/mantis/actions/workflows/typecheck.yml/badge.svg)](https://github.com/uliw/mantis/actions/workflows/typecheck.yml)
[![Documentation Status](https://readthedocs.org/projects/mantis-rtm/badge/?version=latest)](https://mantis-rtm.readthedocs.io/en/latest/?badge=latest)
[![PyPI Version](https://img.shields.io/pypi/v/mantis-rtm.svg)](https://pypi.org/project/mantis-rtm/)
[![Python Versions](https://img.shields.io/badge/python-3.13%20%7C%203.14%20%7C%203.15-blue.svg)](https://pypi.org/project/mantis-rtm/)
[![License: AGPL-3.0-or-later](https://img.shields.io/badge/License-AGPL_3.0--or--later-blue.svg)](https://github.com/uliw/mantis/blob/main/LICENSE)

 Reaction-transport systems involving fast precipitation/dissolution reactions are often extremely stiff. This poses a particular challenge in the geosciences, where the spatial extent of the reaction zone often requires centuries to kiloyears to reach equilibrium, whereas the reaction time scales are often on timescales of seconds. Operator splitting schemes that solve the fast reactions independently of transport are still time-step limited, as larger time steps result in mass imbalances. Here I present a Python framework that uses  a fully coupled monolithic equation system that uses the Newton-Raphson scheme with a Block Thomas Algorithm that can be solved at almost arbitrary large time steps and often within a single iteration.  This allows Mantis to solve stiff non-linear reaction networks within seconds, and in cases where a single iteration is sufficient, milliseconds. Mantis uses 
 an easy-to-read template scheme to define species properties, kinetic constants and reactions including native support for kinetic and equilibrium isotope effects. 


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


## Installation & Documentation 

Please refer to https://mantis-rtm.readthedocs.io/en/latest

## License
This project is licensed under the GPL-3.0 License. See the [LICENSE](LICENSE) file for details.
