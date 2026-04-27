# Example: run a mesh-refinement study for the mixed Dirichlet-Neumann
# diffusion problem with the LDG discretization.

using DiffuSEM
using Trixi

elixir = joinpath(dirname(@__FILE__), "elixir_diffusion_1d_mixed_dirichlet_neumann.jl")

diffusivity = 0.5
forcing_amplitude = 0.4
forcing_frequency = 4.0
dirichlet_mean = 1.0
polydeg = 3
tspan = (0.0, 0.25)
base_initial_refinement_level = 2
iterations = 4
dt_factor = 0.01

println("1D diffusion equation convergence study (mixed Dirichlet-Neumann, LDG)")
println("diffusivity = $(diffusivity), polydeg = $(polydeg), tspan = $(tspan)")
println("forcing_amplitude = $(forcing_amplitude), forcing_frequency = $(forcing_frequency)")
println()

_, errors = convergence_test(@__MODULE__,
                             elixir,
                             iterations;
                             diffusivity = diffusivity,
                             forcing_amplitude = forcing_amplitude,
                             forcing_frequency = forcing_frequency,
                             dirichlet_mean = dirichlet_mean,
                             polydeg = polydeg,
                             tspan = tspan,
                             initial_refinement_level = base_initial_refinement_level,
                             dt_factor = dt_factor,
                             solver_parabolic = ParabolicFormulationLocalDG())

levels = base_initial_refinement_level:(base_initial_refinement_level + iterations - 1)
ndofs = (polydeg + 1) .* (2 .^ collect(levels))
l2_errors = vec(errors[:l2])
linf_errors = vec(errors[:linf])
