using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Celia et al. New Mexico Richards benchmark

problem = HydrologicProblemCeliaNewMexico(tspan = (0.0, 86_400.0))

# Spatial discretization
polydeg = 3
mesh = TreeMesh(problem.domain..., initial_refinement_level = 6, periodicity = false)
solver = DGSEM(; polydeg)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  passive_variables = PassiveVariablesBoundaryFlux1D(),
                                  form = MixedForm())

###############################################################################
# ODE solvers and callbacks

ode = semidiscretize(semi, problem.tspan; jacobian = SparseJacobian())

summary_callback = SummaryCallback()

analysis_interval = 20
analysis_callback = AnalysisCallbackFullState(ode; interval = analysis_interval,
                                              analysis_polydeg = polydeg,
                                              save_analysis = false,
                                              output_directory = "out",
                                              analysis_filename = "analysis.dat",
                                              analysis_errors = Symbol[],
                                              full_state_analysis_integrals =
                                              (water_content, mass_balance))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

# Configure optional spatial adaptivity
amr = true
amr_interval = 10

if amr
    amr_indicator = IndicatorTotalVariation(semi; variable = water_content,
                                            normalize = true)
    amr_controller = ControllerTwoThreshold(semi, amr_indicator;
                                            base_level = 2,
                                            coarsen_threshold = 0.003,
                                            max_level = 10,
                                            refine_threshold = 0.03)
    amr_callback = AMRCallback(semi, amr_controller; interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)
    callbacks = CallbackSet(summary_callback, amr_callback, analysis_callback,
                            alive_callback)
else
    callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)
end

###############################################################################
# run the simulation

run_simulation = true

if run_simulation
    # The mixed form controls error directly in its stored water-content block.
    sol = solve_implicit(ode; dt = 1.0e-2, adaptive = true,
                         reltol = 1.0e-7, abstol = 1.0e-11,
                         saveat = Float64[],
                         error_control_block = evolved_variable_block,
                         error_control_mapping = nothing,
                         isoutofdomain = pressure_head_out_of_domain,
                         callback = callbacks)
end
