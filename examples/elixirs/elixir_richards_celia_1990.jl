using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Celia et al. Richards equation benchmark

problem = HydrologicProblemCelia1990(tspan = (0.0, 360.0))

# Spatial discretization
mesh = TreeMesh(problem.domain..., initial_refinement_level = 5, periodicity = false)
solver = DGSEM(polydeg = 3)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  passive_variables = PassiveVariablesBoundaryFlux1D(),
                                  form = MixedForm())

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, problem.tspan; jacobian = SparseJacobian())

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

if amr
    amr_indicator = IndicatorLoehner(semi, variable = water_content)
    amr_controller = ControllerThreeLevel(semi, amr_indicator; base_level = 1,
                                          med_level = 3, med_threshold = 1.0e-3,
                                          max_level = 6, max_threshold = 1.0e-2)
    amr_callback = AMRCallback(semi, amr_controller; interval = 30,
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
    sol = solve(ode, default_algorithm(semi); dt = 1.0e-2, adaptive = true,
                reltol = 1.0e-7, abstol = 1.0e-11, saveat = Float64[],
                ode_default_options()..., callback = callbacks, maxiters = typemax(Int))
end
