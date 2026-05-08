using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the linear diffusion equation

diffusivity = 0.5
equations = Trixi.LinearDiffusionEquation1D(diffusivity)

initial_refinement_level = 3
polydeg = 3
n_cells_max = 30_000
solver_parabolic = ParabolicFormulationLocalDG()

mesh = TreeMesh((0.0,),
                (1.0,),
                initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max,
                periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)

exact_solution(x, t) = exp(-diffusivity * pi^2 * t) * sinpi(x[1])
initial_condition(x, t, equations) = SVector(exact_solution(x, t))
zero_dirichlet(x, t, equations) = SVector(0.0)

use_boundary_penalty = false
penalty_prefactor = 1.0

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

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.25)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback)

###############################################################################
# run the simulation

dt_factor = 0.01
dt = dt_factor * (1.0 / (2.0^initial_refinement_level))^2
saveat = Float64[]

sol = solve(ode, default_algorithm(semi);
            dt = dt,
            adaptive = false,
            saveat = saveat,
            ode_default_options()..., callback = callbacks,
            maxiters = typemax(Int))
