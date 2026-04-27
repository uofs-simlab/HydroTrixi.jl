using Test
using DiffuSEM
using OrdinaryDiffEq
using SciMLBase
using Trixi

const SMOKE_TESTS = lowercase(get(ENV, "DIFFUSEM_SMOKE_TESTS", "false")) in ("1", "true",
                                                                             "yes", "on")

@testset "1D diffusion equation with Dirichlet BCs" begin
    @test DiffuSEM.examples_dir() == joinpath(pkgdir(DiffuSEM), "examples")

    @test DiffuSEM.default_algorithm() isa OrdinaryDiffEq.Tsit5
    equations = DiffuSEM.Trixi.LinearDiffusionEquation1D(0.1)
    @test equations isa DiffuSEM.Trixi.LinearDiffusionEquation1D

    solver = DGSEM(polydeg = 2, surface_flux = flux_central)
    mesh = TreeMesh((0.0,),
                    (1.0,),
                    initial_refinement_level = 2,
                    n_cells_max = 30_000,
                    periodicity = false)
    initial_condition(x, t, equations) = Trixi.SVector(sinpi(x[1]))

    boundary_conditions = (;
                           x_neg = BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0)),
                           x_pos = BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0)),)

    semi = DiffuSEM.SemidiscretizationParabolic(mesh,
                                                equations,
                                                initial_condition,
                                                solver;
                                                boundary_conditions = boundary_conditions,
                                                solver_parabolic = ParabolicFormulationLocalDG(),)

    @test semi isa DiffuSEM.SemidiscretizationParabolic
    @test semi.solver_parabolic isa Trixi.ParabolicFormulationLocalDG

    ode = semidiscretize(semi, (0.0, 0.01))
    sol = SciMLBase.solve(ode,
                          DiffuSEM.default_algorithm();
                          dt = 1.0e-4,
                          save_everystep = false,)
    @test sol.t[end] ≈ 0.01
    @test all(isfinite, sol.u[end])

    boundary_conditions_mixed = (;
                                 x_neg = BoundaryConditionDirichlet((x, t, equations) -> SVector(1.0)),
                                 x_pos = BoundaryConditionNeumann((x, t, equations) -> SVector(0.0)),)
    semi_mixed = DiffuSEM.SemidiscretizationParabolic(mesh,
                                                      equations,
                                                      initial_condition,
                                                      solver;
                                                      boundary_conditions = boundary_conditions_mixed,
                                                      solver_parabolic = ParabolicFormulationLocalDG(),)
    @test semi_mixed.boundary_conditions.x_neg isa Trixi.BoundaryConditionDirichlet
    @test semi_mixed.boundary_conditions.x_pos isa Trixi.BoundaryConditionNeumann

    ode_mixed = semidiscretize(semi_mixed, (0.0, 0.01))
    sol_mixed = SciMLBase.solve(ode_mixed,
                                DiffuSEM.default_algorithm();
                                dt = 1.0e-4,
                                save_everystep = false,)
    @test sol_mixed.t[end] ≈ 0.01
    @test all(isfinite, sol_mixed.u[end])

    h = 1.0 / (2.0^2)
    penalty_strength = equations.diffusivity * (2 + 1)^2 / h
    boundary_conditions_penalized = (;
                                     x_neg = BoundaryConditionDirichletPenalty((x, t, equations) -> SVector(0.0);
                                                                               penalty = penalty_strength,),
                                     x_pos = BoundaryConditionDirichletPenalty((x, t, equations) -> SVector(0.0);
                                                                               penalty = penalty_strength,),)
    semi_penalized = DiffuSEM.SemidiscretizationParabolic(mesh,
                                                          equations,
                                                          initial_condition,
                                                          solver;
                                                          boundary_conditions = boundary_conditions_penalized,
                                                          solver_parabolic = ParabolicFormulationLocalDG(),)
    ode_penalized = semidiscretize(semi_penalized, (0.0, 0.01))
    sol_penalized = SciMLBase.solve(ode_penalized,
                                    DiffuSEM.default_algorithm();
                                    dt = 1.0e-4,
                                    save_everystep = false,)
    @test sol_penalized.t[end] ≈ 0.01
    @test all(isfinite, sol_penalized.u[end])

    struct ConstantHydraulicConductivity end
    DiffuSEM.hydraulic_conductivity(psi, ::ConstantHydraulicConductivity) = 0.1

    richards_equations = DiffuSEM.RichardsEquation1D(hydraulic_conductivity_model = ConstantHydraulicConductivity())
    constitutive_relation(psi, equations) = psi
    semi_richards = DiffuSEM.SemidiscretizationParabolic(mesh,
                                                         richards_equations,
                                                         initial_condition,
                                                         solver;
                                                         boundary_conditions = boundary_conditions,
                                                         solver_parabolic = ParabolicFormulationLocalDG(),)
    semi_constitutive = DiffuSEM.SemidiscretizationConstitutive(semi_richards,
                                                                constitutive_relation)

    constitutive_initial_state = DiffuSEM.Trixi.compute_coefficients(0.0, semi_constitutive)
    half = length(constitutive_initial_state) ÷ 2
    theta_initial = constitutive_initial_state[1:half]
    psi_initial = constitutive_initial_state[(half + 1):end]
    @test theta_initial ≈ psi_initial

    ode_constitutive = semidiscretize(semi_constitutive, (0.0, 0.01))
    @test size(ode_constitutive.f.mass_matrix) == (2 * half, 2 * half)
    @test all(ode_constitutive.f.mass_matrix[i, i] == 1 for i in 1:half)
    @test all(ode_constitutive.f.mass_matrix[i, i] == 0 for i in (half + 1):(2 * half))

    sol_constitutive = SciMLBase.solve(ode_constitutive,
                                       DiffuSEM.default_algorithm(semi_constitutive);
                                       dt = 1.0e-4,
                                       adaptive = false,
                                       save_everystep = false,)
    @test sol_constitutive.t[end] ≈ 0.01
    @test all(isfinite, sol_constitutive.u[end])

    theta_final = sol_constitutive.u[end][1:half]
    psi_final = sol_constitutive.u[end][(half + 1):end]
    @test theta_final≈psi_final atol=1.0e-9 rtol=1.0e-9

    # With constant conductivity and Θ(ψ) = ψ, the Richards operator recovers the
    # linear diffusion operator in the interior. We check this in a periodic setup
    # so that the gravitational flux offset does not need to be absorbed into
    # boundary data.
    periodic_mesh = TreeMesh((0.0,),
                             (1.0,),
                             initial_refinement_level = 2,
                             n_cells_max = 30_000,
                             periodicity = true)
    periodic_boundary_conditions = Trixi.boundary_condition_periodic
    semi_periodic = DiffuSEM.SemidiscretizationParabolic(periodic_mesh,
                                                         equations,
                                                         initial_condition,
                                                         solver;
                                                         boundary_conditions = periodic_boundary_conditions,
                                                         solver_parabolic = ParabolicFormulationLocalDG(),)
    semi_richards_periodic = DiffuSEM.SemidiscretizationParabolic(periodic_mesh,
                                                                  richards_equations,
                                                                  initial_condition,
                                                                  solver;
                                                                  boundary_conditions = periodic_boundary_conditions,
                                                                  solver_parabolic = ParabolicFormulationLocalDG(),)
    semi_constitutive_periodic = DiffuSEM.SemidiscretizationConstitutive(semi_richards_periodic,
                                                                         constitutive_relation)
    u_diffusion = Trixi.compute_coefficients(0.0, semi_periodic)
    du_diffusion = similar(u_diffusion)
    Trixi.rhs_parabolic!(du_diffusion, u_diffusion, semi_periodic, 0.0)

    u_constitutive = Trixi.compute_coefficients(0.0, semi_constitutive_periodic)
    du_constitutive = similar(u_constitutive)
    ode_constitutive_periodic = semidiscretize(semi_constitutive_periodic, (0.0, 0.01))
    ode_constitutive_periodic.f(du_constitutive, u_constitutive,
                                semi_constitutive_periodic, 0.0)
    periodic_half = length(u_constitutive) ÷ 2
    @test du_constitutive[1:periodic_half]≈du_diffusion atol=1.0e-12 rtol=1.0e-12
    @test du_constitutive[(periodic_half + 1):end]≈zero(du_diffusion) atol=1.0e-12 rtol=1.0e-12

    eoc = DiffuSEM.compute_eoc([1.0, 0.25, 0.0625])
    @test isnan(eoc[1])
    @test eoc[2] ≈ 2.0
    @test eoc[3] ≈ 2.0

    if !SMOKE_TESTS
        mktempdir() do dir
            convergence_plot = joinpath(dir, "convergence.pdf")
            DiffuSEM.plot_convergence_1d([16, 32, 64],
                                         [1.0e-2, 1.0e-3, 1.0e-4],
                                         [2.0e-2, 3.0e-3, 4.0e-4];
                                         output_path = convergence_plot,
                                         triangle_order = 3,)
            @test isfile(convergence_plot)

            solution_plot = joinpath(dir, "solution.pdf")
            DiffuSEM.plot_solution_1d(sol;
                                      output_path = solution_plot,
                                      exact_solution = (x, t) -> exp(-0.1 * pi^2 * t) *
                                                                 sinpi(x[1]),)
            @test isfile(solution_plot)

            constitutive_plot = joinpath(dir, "constitutive_solution.pdf")
            DiffuSEM.plot_solution_1d(sol_constitutive;
                                      output_path = constitutive_plot,
                                      component = 2,)
            @test isfile(constitutive_plot)

            solution_animation = joinpath(dir, "solution.mp4")
            DiffuSEM.animate_solution_1d(sol;
                                         output_path = solution_animation,
                                         framerate = 2,)
            @test isfile(solution_animation)
            @test filesize(solution_animation) > 0

            constitutive_animation = joinpath(dir, "constitutive_solution.mp4")
            DiffuSEM.animate_solution_1d(sol_constitutive;
                                         output_path = constitutive_animation,
                                         component = 1,
                                         framerate = 2,)
            @test isfile(constitutive_animation)
            @test filesize(constitutive_animation) > 0
        end
    end
end
