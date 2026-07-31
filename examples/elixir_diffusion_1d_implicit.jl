using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# implicit semidiscretization of the linear diffusion equation

diffusivity = 0.5
equations = Trixi.LinearDiffusionEquation1D(diffusivity)

initial_refinement_level = 3
polydeg = 3
n_cells_max = 30_000
solver_parabolic = ParabolicFormulationLocalDG()

mesh = TreeMesh((0.0,), (1.0,), initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max, periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)

exact_solution(x, t) = exp(-diffusivity * pi^2 * t) * sinpi(x[1])
initial_condition(x, t, equations) = SVector(exact_solution(x, t))
zero_dirichlet(x, t, equations) = SVector(0.0)
boundary_condition = BoundaryConditionDirichlet(zero_dirichlet)
boundary_conditions = (; x_neg = boundary_condition, x_pos = boundary_condition)

semi_base = SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions = boundary_conditions,
                                        solver_parabolic = solver_parabolic)
semi = SemidiscretizationImplicit(semi_base, TemporalOperatorStandard())

###############################################################################
# ODE solver and callbacks

tspan = (0.0, 0.25)
jacobian = SparseJacobian()
ode = semidiscretize(semi, tspan; jacobian = jacobian)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

###############################################################################
# run the simulation

dt = 1.0e-2
adaptive = true
reltol = 1.0e-9
abstol = 1.0e-11
saveat = Float64[]

sol = solve(ode, default_algorithm(semi); dt = dt, adaptive = adaptive,
            reltol = reltol, abstol = abstol, saveat = saveat,
            ode_default_options()..., callback = callbacks, maxiters = typemax(Int))
