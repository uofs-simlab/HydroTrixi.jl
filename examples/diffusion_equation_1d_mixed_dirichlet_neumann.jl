using DiffuSEM
using SciMLBase
using Trixi

function mixed_dirichlet_neumann_case(;
                                      diffusivity=0.5,
                                      tspan=(0.0, 1.0),
                                      forcing_amplitude=0.4,
                                      forcing_frequency=4.0,
                                      dirichlet_mean=1.0,
                                      parabolic_scheme=ViscousFormulationLocalDG(),
                                      initial_refinement_level=3,
                                      polydeg=3,
                                      n_cells_max=30_000)
    angular_frequency = 2 * pi * forcing_frequency
    mesh = TreeMesh((0.0,), (1.0,),
                    initial_refinement_level=initial_refinement_level,
                    n_cells_max=n_cells_max,
                    periodicity=false)
    solver = DGSEM(polydeg=polydeg, surface_flux=flux_central)
    equations = LinearDiffusionEquation1D(diffusivity)

    complex_wavenumber = sqrt(im * angular_frequency / diffusivity)
    harmonic_shape(x) = cosh(complex_wavenumber * (1 - x)) / cosh(complex_wavenumber)
    exact_solution(x, t) = dirichlet_mean +
                           imag(forcing_amplitude * exp(im * angular_frequency * t) *
                                harmonic_shape(x[1]))

    initial_condition(x, t, equations) = SVector(exact_solution(x, t))
    dirichlet_left(x, t) = dirichlet_mean + forcing_amplitude * sin(angular_frequency * t)
    neumann_right(x, t) = 0.0

    boundary_conditions = (; x_neg=BoundaryConditionDirichlet((x, t, equations) -> SVector(dirichlet_left(x, t))),
                           x_pos=BoundaryConditionNeumann((x, t, equations) -> SVector(neumann_right(x, t))))

    semi = SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                       boundary_conditions=boundary_conditions,
                                       parabolic_scheme=parabolic_scheme)
    ode = semidiscretize(semi, tspan)

    return (; diffusivity, tspan, forcing_amplitude, forcing_frequency, dirichlet_mean,
            angular_frequency, mesh, equations, solver, semi, ode, exact_solution)
end

function solve_mixed_dirichlet_neumann_case(case;
                                            dt=5.0e-5,
                                            adaptive=false,
                                            saveat=nothing,
                                            save_everystep=false,
                                            kwargs...)
    if isnothing(saveat)
        return solve(case.ode, default_algorithm();
                     dt=dt,
                     adaptive=adaptive,
                     save_everystep=save_everystep,
                     kwargs...)
    end

    return solve(case.ode, default_algorithm();
                 dt=dt,
                 adaptive=adaptive,
                 saveat=saveat,
                 kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    case = mixed_dirichlet_neumann_case()
    sol = solve_mixed_dirichlet_neumann_case(case;
                                             dt=5.0e-5,
                                             adaptive=false,
                                             save_everystep=false)
    @show sol.t[end]
end
