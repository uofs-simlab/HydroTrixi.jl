using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Celia et al. Haverkamp Richards benchmark

problem = HydrologicProblemCeliaHaverkamp(tspan = (0.0, 360.0))

# Spatial discretization
mesh = TreeMesh(problem.domain..., initial_refinement_level = 6, periodicity = false)
solver = DGSEM(polydeg = 3)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  passive_variables = PassiveVariablesBoundaryFlux1D(),
                                  form = MixedForm())

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, problem.tspan; jacobian = SparseJacobian())
internalnorm = evolved_variable_norm(semi)
# Control error in water content for both mixed and pressure-head forms.
error_control_variables = water_content

summary_callback = SummaryCallback()

analysis_interval = 20
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     output_directory = "out",
                                     analysis_filename = "analysis.dat",
                                     analysis_errors = Symbol[],
                                     extra_analysis_integrals = (water_content, mass_bias))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

# Configure optional spatial adaptivity
amr = false
amr_interval = 10
amr_base_level = 2

if amr
    amr_indicator = IndicatorTotalVariation(semi; variable = water_content,
                                            normalize = true)
    amr_controller = ControllerTwoThreshold(semi, amr_indicator;
                                            base_level = amr_base_level,
                                            coarsen_threshold = 0.003,
                                            max_level = 10,
                                            refine_threshold = 0.03)
    amr_callback = AMRCallback(semi, amr_controller; interval = amr_interval,
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
                         reltol = 1.0e-7, abstol = 1.0e-11,
                         saveat = Float64[], ode_default_options()...,
                         internalnorm = internalnorm,
                         error_control_variables = error_control_variables,
                         isoutofdomain = pressure_head_out_of_domain,
                         callback = callbacks,
                         maxiters = typemax(Int))
end
