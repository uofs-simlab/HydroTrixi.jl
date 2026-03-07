# DiffuSEM.jl

DiffuSEM.jl is a discontinuous spectral-element framework for solving linear and nonlinear 
parabolic PDEs using the Trixi.jl and SciML ecosystems. The package is scaffolded to 
integrate:

- `Trixi.jl` for semidiscretization setup
- `OrdinaryDiffEq.jl` + `SciMLBase.jl` for time integration

These are required dependencies of `DiffuSEM.jl`. The package architecture is designed to 
support additional linear/nonlinear parabolic PDE models as they are added.

## Package structure

- `src/DiffuSEM.jl`: main module + exports
- `src/equations/`: equation definitions (parabolic PDE models)
- `src/semidiscretization/`: problem and semidiscretization APIs
- `src/time_integration/`: default time-integration algorithm selection
- `src/plotting/solution_plot.jl`: solution plotting routines
- `src/plotting/solution_animation.jl`: solution animation routines
- `src/plotting/convergence_plot.jl`: convergence plotting routines
- `examples/diffusion_equation_1d_mixed_dirichlet_neumann.jl`: reusable mixed-BC case setup + solver helper
- `examples/diffusion_equation_1d_mixed_dirichlet_neumann_plot.jl`: mixed-BC solution plot
- `examples/diffusion_equation_1d_mixed_dirichlet_neumann_animation.jl`: mixed-BC animation
- `examples/diffusion_equation_1d_mixed_dirichlet_neumann_convergence.jl`: mixed-BC convergence study
- `test/runtests.jl`: 1D diffusion-equation tests
