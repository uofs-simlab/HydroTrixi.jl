# HydroTrixi.jl

[![CI](https://github.com/uofs-simlab/HydroTrixi.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/uofs-simlab/HydroTrixi.jl/actions/workflows/ci.yml)
[![docs-dev](https://img.shields.io/badge/docs-dev-blueviolet.svg)](https://tjbmontoya.com/HydroTrixi.jl/dev)

**HydroTrixi.jl** is an adaptive discontinuous spectral-element solver for hydrologic problems. It builds upon the parabolic spatial discretization capabilities in [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) [(Ranocha et al. 2022)](#references) and the time integration methods in [SciML](https://sciml.ai/), adding the following technical features to support the solution of the **Richards equation** in one spatial dimension:

- Implicit temporal operator formulations using algebraic constitutive constraints or capacity functions to support mixed and pressure-head formulations, respectively
- Flexible parabolic boundary conditions allowing penalty-type numerical fluxes that depend on the inner and outer solution and flux values
- Adaptivity in time for stiff systems using high-order embedded Rosenbrock-Wanner methods with sparse Jacobian evaluation
- Adaptivity in space with mass-conservative solution transfer
- Solver flux output method (SFOM) diagnostics for time-integrated boundary fluxes, following [Ireson et al. (2023)](#references)
- Tools to facilitate problem setup and visualization, with examples for standard hydrologic benchmark cases

The mixed formulation advances water content as the conserved variable while evaluating the nonlinear diffusion operator in terms of pressure head, and guarantees a fully discrete water mass balance under spatial as well as temporal adaptivity.

## Installation

If you have not yet installed Julia, please [follow the instructions for your operating system](https://julialang.org/downloads/platform/). HydroTrixi.jl works with Julia v1.10 and newer. We recommend using the latest stable release.

Install and run HydroTrixi.jl from a local clone:

```bash
git clone https://github.com/uofs-simlab/HydroTrixi.jl.git
cd HydroTrixi.jl
julia --project=.
```

Then instantiate dependencies in the Julia REPL:

```julia
julia> using Pkg

julia> Pkg.instantiate()
```

Run the test suite:

```julia
julia> Pkg.test()
```

Worked examples and tutorials are provided in the [documentation](https://tjbmontoya.com/HydroTrixi.jl/dev).

## References

- Ranocha, H., Schlottke-Lakemper, M., Winters, A. R., Faulhaber, E., Chan, J., and Gassner, G. J. (2022). [Adaptive numerical simulations with Trixi.jl: A case study of Julia for scientific computing](https://doi.org/10.21105/jcon.00077). *Proceedings of the JuliaCon Conferences*, 1(1), 77.
- Ireson, A. M., Spiteri, R. J., Clark, M. P., and Mathias, S. A. (2023). [A simple, efficient, mass-conservative approach to solving the Richards equation (openRE, v1.0)](https://doi.org/10.5194/gmd-16-659-2023). *Geoscientific Model Development*, 16, 659-677.

## Acknowledgements

The developers of this package acknowledge funding support from the [Cooperative Institute for Research to Operations in Hydrology (CIROH)](https://ciroh.ua.edu/).

## License

HydroTrixi.jl is released under the [MIT License](https://github.com/uofs-simlab/HydroTrixi.jl/blob/main/LICENSE).
