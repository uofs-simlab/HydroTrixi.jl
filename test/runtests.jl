using Test
using DiffuSEM
using OrdinaryDiffEq
using SciMLBase
using Trixi

@testset "1D diffusion equation with Dirichlet BCs" begin
    @test DiffuSEM.examples_dir() == joinpath(pkgdir(DiffuSEM), "examples")

    @test DiffuSEM.default_algorithm() isa OrdinaryDiffEq.Tsit5
    equations = DiffuSEM.LinearDiffusionEquation1D(0.1)
    @test equations isa DiffuSEM.LinearDiffusionEquation1D

    solver = DGSEM(polydeg=2, surface_flux=flux_central)
    mesh = TreeMesh((0.0,), (1.0,),
                    initial_refinement_level=2,
                    n_cells_max=30_000,
                    periodicity=false)
    initial_condition(x, t, equations) = Trixi.SVector(sinpi(x[1]))

    boundary_conditions = (; x_neg=BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0)),
                           x_pos=BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0)))

    semi = DiffuSEM.SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                                boundary_conditions=boundary_conditions,
                                                parabolic_scheme=ViscousFormulationLocalDG())

    @test semi isa DiffuSEM.SemidiscretizationParabolic
    @test semi.parabolic_scheme isa Trixi.ViscousFormulationLocalDG

    ode = semidiscretize(semi, (0.0, 0.01))
    sol = SciMLBase.solve(ode, DiffuSEM.default_algorithm(); dt=1.0e-4, save_everystep=false)
    @test sol.t[end] ≈ 0.01
    @test all(isfinite, sol.u[end])

    boundary_conditions_mixed = (; x_neg=BoundaryConditionDirichlet((x, t, equations) -> SVector(1.0)),
                                 x_pos=BoundaryConditionNeumann((x, t, equations) -> SVector(0.0)))
    semi_mixed = DiffuSEM.SemidiscretizationParabolic(mesh, equations, initial_condition, solver;
                                                      boundary_conditions=boundary_conditions_mixed,
                                                      parabolic_scheme=ViscousFormulationLocalDG())
    @test semi_mixed.boundary_conditions.x_neg isa Trixi.BoundaryConditionDirichlet
    @test semi_mixed.boundary_conditions.x_pos isa Trixi.BoundaryConditionNeumann

    ode_mixed = semidiscretize(semi_mixed, (0.0, 0.01))
    sol_mixed = SciMLBase.solve(ode_mixed, DiffuSEM.default_algorithm(); dt=1.0e-4,
                                save_everystep=false)
    @test sol_mixed.t[end] ≈ 0.01
    @test all(isfinite, sol_mixed.u[end])

    eoc = DiffuSEM.compute_eoc([1.0, 0.25, 0.0625])
    @test isnan(eoc[1])
    @test eoc[2] ≈ 2.0
    @test eoc[3] ≈ 2.0

    mktempdir() do dir
        convergence_plot = joinpath(dir, "convergence.pdf")
        DiffuSEM.plot_convergence_1d([16, 32, 64], [1.0e-2, 1.0e-3, 1.0e-4],
                                     [2.0e-2, 3.0e-3, 4.0e-4];
                                     output_path=convergence_plot,
                                     triangle_order=3)
        @test isfile(convergence_plot)

        solution_plot = joinpath(dir, "solution.pdf")
        DiffuSEM.plot_solution_1d(sol;
                                  output_path=solution_plot,
                                  exact_solution=(x, t) -> exp(-0.1 * pi^2 * t) * sinpi(x[1]))
        @test isfile(solution_plot)

        solution_animation = joinpath(dir, "solution.mp4")
        DiffuSEM.animate_solution_1d(sol;
                                     output_path=solution_animation,
                                     framerate=2)
        @test isfile(solution_animation)
        @test filesize(solution_animation) > 0
    end
end
