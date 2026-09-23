using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the Richards equation manufactured solution

problem = HydrologicProblemRichardsManufacturedSolution(tspan = (0.0, 120.0))

# Spatial discretization
polydeg = 3
mesh = TreeMesh(problem.domain...; initial_refinement_level = 4, periodicity = false)
solver = DGSEM(; polydeg)

semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  form = MixedForm())

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, problem.tspan; jacobian = SparseJacobian())
algorithm = default_algorithm(ode)

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi; interval = analysis_interval,
                                     analysis_polydeg = polydeg)

callbacks = CallbackSet(SummaryCallback(), analysis_callback,
                        AliveCallback(analysis_interval = analysis_interval))

###############################################################################
# run the simulation

run_simulation = true
solve_options = NamedTuple()

if run_simulation
    sol = solve_implicit(ode, algorithm; dt = 1.0e-2, adaptive = true,
                         reltol = 1.0e-9, abstol = 1.0e-11,
                         saveat = Float64[], callback = callbacks,
                         error_control_block = evolved_variable_block,
                         error_control_mapping = nothing,
                         solve_options...)
    (; sol, analysis_callback)
end
