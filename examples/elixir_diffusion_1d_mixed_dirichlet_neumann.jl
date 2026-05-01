# Example: Trixi-style elixir for the one-dimensional diffusion problem with a
# time-dependent Dirichlet condition at the left boundary and a homogeneous
# Neumann condition at the right boundary.

using HydroTrixi
using SciMLBase
using Trixi

diffusivity = 0.5
tspan = (0.0, 1.0)
forcing_amplitude = 0.4
forcing_frequency = 4.0
dirichlet_mean = 1.0
solver_parabolic = ParabolicFormulationLocalDG()
initial_refinement_level = 3
polydeg = 3
n_cells_max = 30_000
dt_factor = 0.0032
dt = dt_factor * (1.0 / (2.0^initial_refinement_level))^2
adaptive = false
saveat = nothing
save_everystep = false

angular_frequency = 2 * pi * forcing_frequency
mesh = TreeMesh((0.0,),
                (1.0,),
                initial_refinement_level = initial_refinement_level,
                n_cells_max = n_cells_max,
                periodicity = false)
solver = DGSEM(polydeg = polydeg, surface_flux = flux_central)
equations = Trixi.LinearDiffusionEquation1D(diffusivity)

complex_wavenumber = sqrt(im * angular_frequency / diffusivity)
harmonic_shape(x) = cosh(complex_wavenumber * (1 - x)) / cosh(complex_wavenumber)
function exact_solution(x, t)
    dirichlet_mean +
    imag(forcing_amplitude * exp(im * angular_frequency * t) *
         harmonic_shape(x[1]))
end

initial_condition(x, t, equations) = SVector(exact_solution(x, t))
dirichlet_left(x, t) = dirichlet_mean + forcing_amplitude * sin(angular_frequency * t)
neumann_right(x, t) = 0.0

boundary_conditions = (;
                       x_neg = BoundaryConditionDirichlet((x, t, equations) -> SVector(dirichlet_left(x,
                                                                                                      t))),
                       x_pos = BoundaryConditionNeumann((x, t, equations) -> SVector(neumann_right(x,
                                                                                                   t))),)

semi = SemidiscretizationParabolic(mesh,
                                   equations,
                                   initial_condition,
                                   solver;
                                   boundary_conditions = boundary_conditions,
                                   solver_parabolic = solver_parabolic)
analysis_callback = AnalysisCallback(semi, interval = typemax(Int))
ode = semidiscretize(semi, tspan)

if isnothing(saveat)
    sol = solve(ode,
                default_algorithm(semi);
                dt = dt,
                adaptive = adaptive,
                save_everystep = save_everystep)
else
    sol = solve(ode,
                default_algorithm(semi);
                dt = dt,
                adaptive = adaptive,
                saveat = saveat)
end
