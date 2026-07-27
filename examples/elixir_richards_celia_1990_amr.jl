using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Celia et al. Richards equation benchmark

final_time = 360.0
problem = HydrologicProblemCelia1990(tspan = (0.0, final_time))

coordinates_min, coordinates_max = problem.domain
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 30_000,
                periodicity = false)
solver = DGSEM(polydeg = 3, surface_flux = flux_central)
passive_variables = PassiveVariablesBoundaryFlux1D()
form = MixedForm()

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  passive_variables = passive_variables,
                                  form = form)

###############################################################################
# ODE solvers, callbacks etc.

tspan = problem.tspan
sparse_jacobian = false
ode = semidiscretize(semi, tspan; sparse_jacobian = sparse_jacobian)

summary_callback = SummaryCallback()

analysis_interval = 10
extra_analysis_integrals = (HydroTrixi.water_content,
                            HydroTrixi.mass_bias)
save_analysis = false
output_directory = "out"
analysis_filename = "analysis.dat"
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = save_analysis,
                                     output_directory = output_directory,
                                     analysis_filename = analysis_filename,
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals = extra_analysis_integrals)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

amr_indicator = IndicatorLoehner(semi, variable = HydroTrixi.water_content)
amr_interval = 30
max_level = 6
amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                      base_level = 1,
                                      med_level = 3, med_threshold = 1e-3,
                                      max_level = max_level, max_threshold = 1e-2)
amr_callback = AMRCallback(semi, amr_controller;
                           interval = amr_interval,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        amr_callback)

###############################################################################
# run the simulation

dt = 1.0e-2
adaptive = true
reltol = 1.0e-7
abstol = 1.0e-11
saveat = Float64[]
run_simulation = true

if run_simulation
    sol = solve(ode, default_algorithm(semi);
                dt = dt,
                adaptive = adaptive,
                reltol = reltol,
                abstol = abstol,
                saveat = saveat,
                ode_default_options()..., callback = callbacks)
end
