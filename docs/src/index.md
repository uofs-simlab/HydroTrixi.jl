```@meta
EditURL = "https://github.com/uofs-simlab/HydroTrixi.jl/blob/main/docs/src/index.md"
```

# HydroTrixi.jl

**HydroTrixi.jl** is an adaptive discontinuous spectral-element solver for hydrologic
problems. It builds upon the parabolic spatial discretization capabilities in
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl) and the time integration
methods in the [SciML ecosystem](https://sciml.ai/), with support for one-dimensional
Richards and diffusion problems.

## Tutorials

- [Celia *et al.* (1990) infiltration problem](tutorials/celia_1990.md)
- [Sparse finite-difference Jacobian evaluation](tutorials/richards_jacobian_sparsity.md)

## Installation

HydroTrixi.jl supports Julia v1.10 and newer. To work with a local checkout:

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

## Acknowledgements

The developers of this package acknowledge funding support from the
[Cooperative Institute for Research to Operations in Hydrology
(CIROH)](https://ciroh.ua.edu/).
