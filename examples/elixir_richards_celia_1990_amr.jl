using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Celia et al. Richards equation benchmark

problem = HydrologicProblemCelia1990()

coordinates_min, coordinates_max = problem.domain
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 30_000,
                periodicity = false)
solver = DGSEM(polydeg = 3, surface_flux = flux_central)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG())

###############################################################################
# ODE solvers, callbacks etc.

tspan = problem.tspan
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals = (HydroTrixi.water_content,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

amr_indicator = IndicatorLoehner(semi, variable = HydroTrixi.water_content)
amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                      base_level = 1,
                                      med_level = 3, med_threshold = 1e-3,
                                      max_level = 7, max_threshold = 1e-2)
amr_callback = AMRCallback(semi, amr_controller;
                           interval = 10,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        amr_callback)

###############################################################################
# run the simulation

dt = 1.0e-2
adaptive = true
saveat = Float64[]
run_simulation = true

if run_simulation
    sol = solve(ode, default_algorithm(semi);
                dt = dt,
                adaptive = adaptive,
                saveat = saveat,
                ode_default_options()..., callback = callbacks)
end
