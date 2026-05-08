using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Celia et al. Richards equation benchmark

problem = HydrologicProblemCelia1990()

initial_refinement_level = 5
polydeg = 3
n_cells_max = 30_000
solver_parabolic = ParabolicFormulationLocalDG()

coordinates_min, coordinates_max = problem.domain
mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max,
                periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = solver_parabolic)

###############################################################################
# ODE solvers, callbacks etc.

tspan = problem.tspan
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 20
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals = (HydroTrixi.water_content,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback)

###############################################################################
# run the simulation

dt = 1.0e-2
adaptive = true
saveat = Float64[]

sol = solve(ode, default_algorithm(semi);
            dt = dt,
            adaptive = adaptive,
            saveat = saveat,
            ode_default_options()..., callback = callbacks,
            maxiters = typemax(Int))
