# Example: Trixi-style elixir for the one-dimensional diffusion problem with
# homogeneous Dirichlet conditions at both boundaries.

using DiffuSEM
using SciMLBase
using Trixi

diffusivity = 0.5
tspan = (0.0, 0.25)
solver_parabolic = ParabolicFormulationLocalDG()
initial_refinement_level = 3
polydeg = 3
n_cells_max = 30_000
dt_factor = 0.01
dt = dt_factor * (1.0 / (2.0^initial_refinement_level))^2
adaptive = false
saveat = nothing
save_everystep = false
use_boundary_penalty = false
penalty_prefactor = 1.0

mesh = TreeMesh((0.0,),
                (1.0,),
                initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max,
                periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)
equations = Trixi.LinearDiffusionEquation1D(diffusivity)

exact_solution(x, t) = exp(-diffusivity * pi^2 * t) * sinpi(x[1])
initial_condition(x, t, equations) = SVector(exact_solution(x, t))
zero_dirichlet(x, t, equations) = SVector(0.0)

if use_boundary_penalty
    h = 1.0 / (2.0^initial_refinement_level)
    penalty = penalty_prefactor * diffusivity * (polydeg + 1)^2 / h
    boundary_conditions = (;
                           x_neg = BoundaryConditionDirichletPenalty(zero_dirichlet;
                                                                     penalty = penalty),
                           x_pos = BoundaryConditionDirichletPenalty(zero_dirichlet;
                                                                     penalty = penalty))
else
    boundary_conditions = (;
                           x_neg = BoundaryConditionDirichlet(zero_dirichlet),
                           x_pos = BoundaryConditionDirichlet(zero_dirichlet))
end

semi = SemidiscretizationParabolic(mesh,
                                   equations,
                                   initial_condition,
                                   solver;
                                   boundary_conditions = boundary_conditions,
                                   solver_parabolic = solver_parabolic)
analysis_callback = AnalysisCallback(semi, interval = typemax(Int))
ode = semidiscretize(semi, tspan)

if isnothing(saveat)
    sol = solve(ode,
                default_algorithm();
                dt = dt,
                adaptive = adaptive,
                save_everystep = save_everystep,
                maxiters = typemax(Int))
else
    sol = solve(ode,
                default_algorithm();
                dt = dt,
                adaptive = adaptive,
                saveat = saveat,
                maxiters = typemax(Int))
end
