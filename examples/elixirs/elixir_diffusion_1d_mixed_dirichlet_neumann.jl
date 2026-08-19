using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the linear diffusion equation

# Set up a standard linear diffusion equation with constant diffusivity
diffusivity = 0.5
equations = Trixi.LinearDiffusionEquation1D(diffusivity)

# Spatial discretization
mesh = TreeMesh((0.0,), (1.0,), initial_refinement_level = 3, periodicity = false)
solver = DGSEM(polydeg = 3)

# Define the exact solution and initial condition
forcing_amplitude = 0.4
angular_frequency = 2 * pi * 4.0
complex_wavenumber = sqrt(im * angular_frequency / diffusivity)
function exact_solution(x, t)
    dirichlet_mean + imag(forcing_amplitude * exp(im * angular_frequency * t) *
         cosh(complex_wavenumber * (1 - x[1])) / cosh(complex_wavenumber))
end
initial_condition(x, t, equations) = SVector(exact_solution(x, t))

# Define the boundary conditions
dirichlet_mean = 1.0
function dirichlet_left_boundary(x, t, equations)
    SVector(dirichlet_mean + forcing_amplitude * sin(angular_frequency * t))
end
neumann_right_boundary(x, t, equations) = SVector(0.0)
boundary_conditions = (; x_neg = BoundaryConditionDirichlet(dirichlet_left_boundary),
                       x_pos = BoundaryConditionNeumann(neumann_right_boundary))

# Construct an identity-mass-matrix representation for the implicit solver
semi_base = SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions = boundary_conditions,
                                        solver_parabolic = ParabolicFormulationLocalDG())
semi = SemidiscretizationImplicit(semi_base, TemporalOperatorStandard())

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

###############################################################################
# run the simulation

sol = solve(ode, default_algorithm(ode);
            dt = 0.0032 / length(Trixi.leaf_cells(mesh.tree))^2, adaptive = false,
            saveat = Float64[], ode_default_options()..., callback = callbacks)
