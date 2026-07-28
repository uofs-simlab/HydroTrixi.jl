```@meta
EditURL = "https://github.com/uofs-simlab/HydroTrixi.jl/blob/main/docs/src/index.md"
```

# HydroTrixi.jl

**HydroTrixi.jl** is an adaptive discontinuous spectral-element solver for hydrologic
problems. It builds upon the parabolic spatial discretization capabilities in
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl) and the time integration
methods in the [SciML ecosystem](https://sciml.ai/), with support for one-dimensional
Richards and diffusion problems.

For the Richards equation, HydroTrixi.jl uses an arbitrary-order local discontinuous
Galerkin spectral-element discretization with collocated Legendre-Gauss-Lobatto
quadrature. The mixed formulation evolves water content directly while enforcing its
constitutive relation with pressure head as an algebraic constraint. The pressure-head
formulation instead advances pressure head through the nonlinear capacity function. Both
forms satisfy the same semi-discrete water mass balance, while direct evolution of water
content enables the mixed formulation to retain that linear balance under time
integration and conservative mesh transfer. Adaptive Rosenbrock-Wanner methods provide
temporal adaptivity, and adaptive mesh refinement can resolve sharp wetting fronts.

## Tutorials

- [Celia *et al.* (1990) infiltration problem](tutorials/celia_1990.md)
- [Sparse finite-difference Jacobian evaluation for the mixed
  formulation](tutorials/richards_jacobian_sparsity.md)

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
