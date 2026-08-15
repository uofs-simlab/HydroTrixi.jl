using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the linear diffusion equation

diffusivity = 0.5
equations = Trixi.LinearDiffusionEquation1D(diffusivity)

# Spatial discretization
mesh = TreeMesh((0.0,), (1.0,), initial_refinement_level = 3, periodicity = false)
solver = DGSEM(polydeg = 3)

exact_solution(x, t) = exp(-diffusivity * pi^2 * t) * sinpi(x[1])
initial_condition(x, t, equations) = SVector(exact_solution(x, t))
zero_dirichlet(x, t, equations) = SVector(0.0)

# Set this factor to zero to recover the unpenalized divergence flux
penalty_factor = 0
boundary_condition = BoundaryConditionDirichletPenalty(zero_dirichlet; penalty_factor)
boundary_conditions = (; x_neg = boundary_condition, x_pos = boundary_condition)

# Construct an identity-mass-matrix representation for explicit and implicit solvers
semi_base = SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions = boundary_conditions,
                                        solver_parabolic = ParabolicFormulationLocalDG())
semi = SemidiscretizationImplicit(semi_base, TemporalOperatorStandard())

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.25)
ode = semidiscretize(semi, tspan; jacobian = DenseJacobian())
algorithm = default_algorithm(semi)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

###############################################################################
# run the simulation

sol = solve(ode, algorithm;
            dt = 0.01 / length(Trixi.leaf_cells(mesh.tree))^2, adaptive = false,
            reltol = 1.0e-9, abstol = 1.0e-11, saveat = Float64[], ode_default_options()...,
            callback = callbacks, maxiters = typemax(Int))
