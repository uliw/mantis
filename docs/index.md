# MANTIS Documentation

**MANTIS** (*A Monolithic Assembly of Nonlinear Transport Interactions Solver*) is a high-performance Python package for reactive transport modeling (RTM) in aquatic and sediment geochemistry.

## Features
- **Monolithic Block Newton Solver**: Coupled implicit solution of advective-diffusive-reactive transport and complex reaction networks.
- **Accurate Speciation & Thermodynamics**: Built-in seawater and porewater carbonate equilibria.
- **Isotope Fractionation**: Automated parsing and solving of multi-isotope systems (e.g. $\delta^{34}\text{S}$, $\delta^{13}\text{C}$).
- **Numba Acceleration**: High-speed vector norms, residual evaluations, and step-control kernels.

```{toctree}
:maxdepth: 2
:caption: User Guide & Tutorials

Tutorial_1
Tutorial_2
```

```{toctree}
:maxdepth: 2
:caption: API Reference

api/index
```

## Indices and Tables
* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`
