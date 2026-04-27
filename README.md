# DiffuSEM.jl

[![CI](https://github.com/tristanmontoya/DiffuSEM.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/tristanmontoya/DiffuSEM.jl/actions/workflows/ci.yml)

DiffuSEM.jl is a discontinuous spectral-element framework for solving linear and nonlinear 
parabolic PDEs using the Trixi.jl and SciML ecosystems. The package architecture is designed to 
support additional linear/nonlinear parabolic PDE models as they are added.

## Installation

If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). DiffuSEM.jl works
with Julia v1.10 and newer. We recommend using the latest stable release.

Install and run DiffuSEM.jl from a local clone:

```bash
git clone https://github.com/tristanmontoya/DiffuSEM.jl.git
cd DiffuSEM.jl
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

## Example usage

### Linear diffusion with mixed Dirichlet-Neumann boundary conditions

Run the simulation using `trixi_include`:

```julia
using DiffuSEM
using Trixi

elixir = joinpath(DiffuSEM.examples_dir(),
                  "elixir_diffusion_1d_mixed_dirichlet_neumann.jl")
tspan = (0.0, 1.0)
trixi_include(elixir;
              tspan = tspan,
              saveat = range(first(tspan), last(tspan); length = 121))
```

Plot the final-time solution profile at `t = 1.0` using the provided
visualization utilities based on `CairoMakie`:

```julia
plot_solution_1d(sol;
                 exact_solution = exact_solution,
                 title = "t = 1.0",
                 output_path = joinpath("plots",
                                        "diffusion_equation_1d_mixed_dirichlet_neumann_solution.png"))
```

![Mixed BC solution profile](assets/images/diffusion_equation_1d_mixed_dirichlet_neumann_solution.png)

Generate an animation from the same solution:
```julia
animate_solution_1d(sol;
                    exact_solution = exact_solution,
                    output_path = joinpath("plots",
                                           "diffusion_equation_1d_mixed_dirichlet_neumann_solution_long.gif"),
                    ylims = (0.5, 1.5),
                    framerate = 30)
```

![Mixed BC solution animation](assets/images/diffusion_equation_1d_mixed_dirichlet_neumann_solution_long.gif)

Run the BR1-vs-LDG convergence study:

```julia
include(joinpath(DiffuSEM.examples_dir(),
                 "convergence_diffusion_1d_mixed_dirichlet_neumann_br1_vs_ldg.jl"))
```

Generate the convergence plot from the computed study data:

```julia
include(joinpath(DiffuSEM.examples_dir(),
                 "plot_convergence_diffusion_1d_mixed_dirichlet_neumann_br1_vs_ldg.jl"))
```

![BR1 vs LDG convergence](assets/images/diffusion_equation_1d_mixed_dirichlet_neumann_convergence_br1_vs_ldg.png)
