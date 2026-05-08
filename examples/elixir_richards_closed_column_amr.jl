using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the closed-column Richards equation

final_time = 86400.0
tspan = (0.0, final_time)

# Use standard loam van Genuchten-Mualem parameters in SI units
soil_model = VanGenuchten(saturated_hydraulic_conductivity = 24.96e-2 / 86_400,
                          alpha = 3.6,
                          n = 1.56,
                          theta_s = 0.43,
                          theta_r = 0.078)

problem = HydrologicProblemRichardsClosedColumn(tspan = tspan,
                                                soil_model = soil_model)

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

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10
extra_analysis_integrals = (HydroTrixi.water_content,
                            HydroTrixi.water_content_timederivative)
analysis_callback = AnalysisCallback(semi,
                                     interval = analysis_interval,
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals = extra_analysis_integrals)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

amr_indicator = IndicatorLoehner(semi, variable = HydroTrixi.water_content)
amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                      base_level = 1,
                                      med_level = 3, med_threshold = 1.0e-3,
                                      max_level = 7, max_threshold = 1.0e-2)
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
                ode_default_options()..., callback = callbacks,
                maxiters = typemax(Int))
end
