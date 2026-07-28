using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Richards equation manufactured solution

final_time = 120.0
problem = HydrologicProblemRichardsManufacturedSolution(tspan = (0.0, final_time))

initial_refinement_level = 4
polydeg = 3
n_cells_max = 30_000
solver_parabolic = ParabolicFormulationLocalDG()
form = MixedForm()

coordinates_min, coordinates_max = problem.domain
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max, periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = solver_parabolic, form = form)

###############################################################################
# ODE solvers, callbacks etc.

tspan = problem.tspan
jacobian = DefaultJacobian()
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
run_simulation = true

if run_simulation
    sol = solve(ode, default_algorithm(semi, jacobian); dt = dt, adaptive = adaptive,
                reltol = reltol, abstol = abstol, saveat = saveat,
                ode_default_options()..., callback = callbacks, maxiters = typemax(Int))
end
