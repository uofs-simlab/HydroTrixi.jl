using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Richards equation manufactured solution

problem = HydrologicProblemRichardsManufacturedSolution(tspan = (0.0, 120.0))

# Spatial discretization
mesh = TreeMesh(problem.domain..., initial_refinement_level = 4, periodicity = false)
solver = DGSEM(polydeg = 3)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  form = MixedForm())

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, problem.tspan; jacobian = SparseJacobian())
algorithm = default_algorithm(semi)
internalnorm = state_variable_norm(semi)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

###############################################################################
# run the simulation

run_simulation = true

if run_simulation
    sol = solve(ode, algorithm; dt = 1.0e-2, adaptive = true,
                reltol = 1.0e-9, abstol = 1.0e-11, saveat = Float64[],
                ode_default_options()..., internalnorm = internalnorm,
                callback = callbacks, maxiters = typemax(Int))
end
