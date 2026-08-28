using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the closed-column Richards equation

# Use standard loam van Genuchten-Mualem parameters in SI units
soil_model = VanGenuchten(saturated_hydraulic_conductivity = 24.96e-2 / 86_400, alpha = 3.6,
                          n = 1.56, theta_s = 0.43, theta_r = 0.078)
problem = HydrologicProblemRichardsClosedColumn(tspan = (0.0, 360.0),
                                                soil_model = soil_model)

# Spatial discretization
mesh = TreeMesh(problem.domain..., initial_refinement_level = 5, periodicity = false)
solver = DGSEM(polydeg = 3)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  form = MixedForm())

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, problem.tspan; jacobian = SparseJacobian())
internalnorm = evolved_variable_norm(semi)

summary_callback = SummaryCallback()

analysis_interval = 10
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals =
                                     (water_content, water_content_timederivative))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

# Configure optional spatial adaptivity
amr = false

if amr
    amr_indicator = IndicatorTotalVariation(semi; variable = effective_saturation)
    amr_controller = ControllerTwoThreshold(semi, amr_indicator; base_level = 1,
                                            coarsen_threshold = 0.03,
                                            max_level = 7, refine_threshold = 0.09)
    amr_callback = AMRCallback(semi, amr_controller; interval = 20,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)
    callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                            amr_callback)
else
    callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)
end

###############################################################################
# run the simulation

run_simulation = true

if run_simulation
    sol = solve_implicit(ode; dt = 1.0e-2, adaptive = true,
                         saveat = Float64[], ode_default_options()...,
                         internalnorm = internalnorm, callback = callbacks,
                         maxiters = typemax(Int))
end
