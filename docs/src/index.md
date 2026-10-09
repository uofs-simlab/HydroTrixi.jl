```@meta
EditURL = "https://github.com/uofs-simlab/HydroTrixi.jl/blob/main/docs/src/index.md"
```

# HydroTrixi.jl

**HydroTrixi.jl** is an adaptive discontinuous spectral-element solver for hydrologic problems. It builds upon the parabolic spatial discretization capabilities in [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) and the time integration methods in the [SciML ecosystem](https://sciml.ai/), with support for the one-dimensional Richards equation.

HydroTrixi.jl uses an arbitrary-order local discontinuous Galerkin spectral-element discretization with collocated Legendre-Gauss-Lobatto quadrature. The mixed formulation evolves water content directly while enforcing its constitutive relation with pressure head as an algebraic constraint. The pressure-head formulation instead advances pressure head through the nonlinear capacity function. For the source-free Richards equation, both formulations satisfy the same semi-discrete water balance. The mixed formulation with Rosenbrock-Wanner integration and water-content transfer also preserves the fully discrete balance in exact arithmetic, with numerical calculations subject to roundoff. The pressure-head formulation has no corresponding fully discrete conservation guarantee. Adaptive Rosenbrock-Wanner methods provide temporal adaptivity, and adaptive mesh refinement can resolve sharp wetting fronts.

## Tutorials

- [Haverkamp infiltration problem](tutorials/celia_haverkamp.md)
- [Sparse Jacobian evaluation](tutorials/richards_jacobian_sparsity.md)

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

The developers of this package acknowledge funding support from the [Cooperative Institute for Research to Operations in Hydrology (CIROH)](https://ciroh.ua.edu/).
