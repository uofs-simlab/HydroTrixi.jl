# HydroTrixi.jl

[![CI](https://github.com/uofs-simlab/HydroTrixi.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/uofs-simlab/HydroTrixi.jl/actions/workflows/ci.yml)
[![docs-dev](https://img.shields.io/badge/docs-dev-blueviolet.svg)](https://tjbmontoya.com/HydroTrixi.jl/dev)

**HydroTrixi.jl** is an adaptive discontinuous spectral-element solver for hydrologic problems. It builds upon the parabolic spatial discretization capabilities in
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl) and the time integration methods in [SciML ecosystem](https://sciml.ai/), adding the following technical features to support the solution of the **Richards equation** in one spatial dimension:
- Mixed formulations that can advance different variables in time from those used in the spatial operator, with the constitutive relation imposed as a constraint as part of a differential algebraic equation
- Flexible parabolic boundary conditions allowing penalty-type numerical fluxes that depend on the inner and outer solution and flux values
- Tools to facilitate problem setup and visualization, with examples for standard hydrologic bencmark cases

## Installation

If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). HydroTrixi.jl works
with Julia v1.10 and newer. We recommend using the latest stable release.

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

Optional sanity check:

```julia
julia> Pkg.test()
```

Visualization is optional. Loading `CairoMakie` and `LaTeXStrings` activates
the plotting extension.

Worked examples and tutorials are provided in the
[documentation](https://tjbmontoya.com/HydroTrixi.jl/dev).

## Acknowledgements
The developers of this package acknowledge funding support from the [Cooperative Institute for Research to Operations in Hydrology (CIROH)](https://ciroh.ua.edu/).
