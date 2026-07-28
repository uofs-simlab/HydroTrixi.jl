using HydroTrixi
using SciMLBase
using Trixi

###############################################################################
# semidiscretization of the linear diffusion equation

diffusivity = 0.5
forcing_amplitude = 0.4
forcing_frequency = 4.0
dirichlet_mean = 1.0
equations = Trixi.LinearDiffusionEquation1D(diffusivity)

initial_refinement_level = 3
polydeg = 3
n_cells_max = 30_000
solver_parabolic = ParabolicFormulationLocalDG()

angular_frequency = 2 * pi * forcing_frequency
mesh = TreeMesh((0.0,), (1.0,), initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max, periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)

complex_wavenumber = sqrt(im * angular_frequency / diffusivity)
function exact_solution(x, t)
    dirichlet_mean + imag(forcing_amplitude * exp(im * angular_frequency * t) *
         cosh(complex_wavenumber * (1 - x[1])) / cosh(complex_wavenumber))
end

initial_condition(x, t, equations) = SVector(exact_solution(x, t))
function dirichlet_left_boundary(x, t, equations)
    SVector(dirichlet_mean + forcing_amplitude * sin(angular_frequency * t))
end
neumann_right_boundary(x, t, equations) = SVector(0.0)

boundary_conditions = (; x_neg = BoundaryConditionDirichlet(dirichlet_left_boundary),
                       x_pos = BoundaryConditionNeumann(neumann_right_boundary))

semi = SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                   boundary_conditions = boundary_conditions,
                                   solver_parabolic = solver_parabolic)

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

dt_factor = 0.0032
dt = dt_factor * (1.0 / (2.0^initial_refinement_level))^2
saveat = Float64[]

sol = solve(ode, default_algorithm(semi); dt = dt, adaptive = false, saveat = saveat,
            ode_default_options()..., callback = callbacks)
