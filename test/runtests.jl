using Test
using HydroTrixi
import OrdinaryDiffEqRosenbrock
import SciMLBase
using SciMLBase: DiscreteCallback, solve, successful_retcode
using Trixi: trixi_include
import Trixi
using TrixiTest

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")

macro test_trixi_include(args...)
    esc(Expr(:macrocall, Symbol("@test_trixi_include_base"), __source__, args...))
end

@testset "sparse-AD Jacobian" begin
    for form in (PressureHeadForm(), MixedForm())
        problem = HydrologicProblemRichardsManufacturedSolution()
        mesh = Trixi.TreeMesh(problem.domain...; initial_refinement_level = 1,
                              periodicity = false)
        semi = SemidiscretizationImplicit(mesh, problem, Trixi.DGSEM(polydeg = 3);
                                          solver_parabolic =
                                          Trixi.ParabolicFormulationLocalDG(),
                                          passive_variables =
                                          PassiveVariablesBoundaryFlux1D(), form)
        sparse_ode = Trixi.semidiscretize(semi, (0.0, 1.0e-3);
                                          jacobian = SparseJacobian())
        dense_ode = Trixi.semidiscretize(semi, (0.0, 1.0e-3);
                                         jacobian = DenseJacobian())
        sparse_integrator = SciMLBase.init(sparse_ode, default_algorithm(sparse_ode);
                                           dt = 1.0e-3, adaptive = false)
        finite_difference = OrdinaryDiffEqRosenbrock.AutoFiniteDiff()
        dense_algorithm = default_algorithm(dense_ode; autodiff = finite_difference)
        dense_integrator = SciMLBase.init(dense_ode, dense_algorithm;
                                          dt = 1.0e-3, adaptive = false)

        SciMLBase.step!(sparse_integrator)
        SciMLBase.step!(dense_integrator)

        @test Matrix(sparse_integrator.cache.J)≈dense_integrator.cache.J rtol=1.0e-6

        # Automatic differentiation rebuilds both spatial caches with dual storage.
        dual_caches = sparse_ode.f.f.dual_semidiscretizations.bufs
        @test !isempty(dual_caches)
        for dual_semi in values(dual_caches)
            @test eltype(dual_semi.cache_parabolic.parabolic_boundaries.flux_values) ==
                  eltype(dual_semi.semi_base.cache.elements)
            @test dual_semi.cache_parabolic.parabolic_container ===
                  dual_semi.semi_base.cache_parabolic.parabolic_container
        end
    end
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl dense Jacobian" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        l2=[7.161872258671082e-5], linf=[0.00035212174570417587])

    @test semi.semi_base === semi_base
    @test semi.cache_parabolic isa HydroTrixi.CacheParabolic1D
    @test semi.cache_parabolic.parabolic_container ===
          semi_base.cache_parabolic.parabolic_container
    @test !hasproperty(semi_base.cache_parabolic, :parabolic_boundaries)
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl sparse Jacobian" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        algorithm=default_algorithm(ode), jacobian=SparseJacobian(),
                        dt=1.0e-2, adaptive=true,
                        reltol=1.0e-9, abstol=1.0e-11,
                        l2=[7.161872258685223e-5], linf=[0.00035212174564494547])
end

@trixi_testset "elixir_diffusion_1d_mixed_dirichlet_neumann.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_mixed_dirichlet_neumann.jl"),
                        dt=1.0e-3,
                        l2=[4.1106696084102e-5], linf=[0.00022679747809317696])
end

@trixi_testset "elixir_diffusion_1d_mixed_dirichlet_neumann.jl native Trixi" begin
    import OrdinaryDiffEqRosenbrock

    # Native spatial caches retain upstream dispatch after HydroTrixi is loaded.
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_mixed_dirichlet_neumann.jl"),
                        semi=semi_base,
                        algorithm=OrdinaryDiffEqRosenbrock.Rodas5P(
                            autodiff=OrdinaryDiffEqRosenbrock.AutoFiniteDiff()),
                        dt=1.0e-3,
                        l2=[4.1106696084102e-5], linf=[0.00022679747809317696])

    @test SciMLBase.successful_retcode(sol)
    @test !hasproperty(semi.cache_parabolic, :parabolic_boundaries)
    u_ode = last(sol.u)
    GC.@preserve u_ode begin
        u = Trixi.wrap_array(u_ode, semi)
        arguments = (u, u, last(sol.t), semi.mesh, semi.equations,
                     semi.boundary_conditions, semi.source_terms, semi.solver,
                     semi.solver_parabolic, semi.cache, semi.cache_parabolic)
        @test which(Trixi.rhs_parabolic!, typeof(arguments)).module === Trixi
    end

    # Non-parabolic discretizations use the native spatial fallback.
    equations_hyperbolic = Trixi.LinearScalarAdvectionEquation1D(1.0)
    semi_hyperbolic = Trixi.SemidiscretizationHyperbolic(mesh, equations_hyperbolic,
                                                        initial_condition, solver;
                                                        boundary_conditions =
                                                        Trixi.BoundaryConditionDirichlet(initial_condition))
    implicit_hyperbolic = SemidiscretizationImplicit(semi_hyperbolic,
                                                     TemporalOperatorStandard())
    @test implicit_hyperbolic.semi_base === semi_hyperbolic
    @test isnothing(implicit_hyperbolic.cache_parabolic)
    state = Trixi.compute_coefficients(0.0, implicit_hyperbolic)
    native_rhs = similar(state)
    implicit_rhs = similar(state)
    Trixi.default_rhs(semi_hyperbolic)(native_rhs, state, semi_hyperbolic, 0.0)
    Trixi.default_rhs(implicit_hyperbolic)(implicit_rhs, state, implicit_hyperbolic, 0.0)
    @test implicit_rhs == native_rhs
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl native 2D AMR" begin
    # The exact profile is independent of y, so it also solves the 2D diffusion equation.
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        mesh=Trixi.TreeMesh((0.0, 0.0), (1.0, 1.0);
                                            initial_refinement_level=1, periodicity=false),
                        equations=Trixi.LinearDiffusionEquation2D(diffusivity),
                        boundary_conditions=Trixi.BoundaryConditionDirichlet(
                            initial_condition),
                        tspan=(0.0, 0.01), dt=0.001,
                        l2=[0.0013200742189360601], linf=[0.004697866874314283])

    @test SciMLBase.successful_retcode(sol)
    @test isnothing(semi.cache_parabolic)
    state = copy(last(sol.u))
    initial_elements = Trixi.nelements(solver, semi_base.cache)
    request = Ref(1)
    controller = (u, mesh, equations, dg, cache; kwargs...) ->
                 fill(request[], Trixi.nelements(dg, cache))
    amr_callback = Trixi.AMRCallback(semi, controller; interval = 1,
                                    adapt_initial_condition = false)

    # Refine every cell, then coarsen complete sibling groups back to the original mesh.
    for (direction, element_factor) in ((1, 4), (-1, 1))
        request[] = direction
        @test amr_callback.affect!(state, semi, last(sol.t), 0)
        n_elements = Trixi.nelements(solver, semi_base.cache)
        @test n_elements == element_factor * initial_elements
        @test size(semi_base.cache_parabolic.parabolic_container.u_transformed, 4) ==
              n_elements
        du = similar(state)
        Trixi.default_rhs(semi)(du, state, semi, last(sol.t))
        @test all(isfinite, du)
    end
    @test state≈last(sol.u) rtol=1.0e-12
end

@trixi_testset "elixir_richards_celia_haverkamp.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        amr=true, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        l2=[0.05128549457329617], linf=[0.41441386613429104])

    @test SciMLBase.successful_retcode(sol)
    @test result isa ImplicitSolveResult
    @test result.sol === sol
    @test isnothing(result.mesh_history)
    @test sol.stats.naccept == 83
    @test sol.stats.nreject > 0
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    levels = Trixi.current_element_levels(mesh, dg, cache)
    @test extrema(levels) == (1, 4)
    @test length(sol.u[end]) > length(ode.u0)

    # Normalization remains finite for constant and nearly constant states.
    u = fill(-0.5, size(cache.elements.node_coordinates))
    @test all(iszero, amr_indicator(u, mesh, equations, dg, cache))
    u[1] += eps()
    @test all(isfinite, amr_indicator(u, mesh, equations, dg, cache))
end

@trixi_testset "elixir_richards_celia_haverkamp.jl accepted-step history" begin
    analysis_directory = mktempdir()
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        amr=true, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        analysis_interval=1, save_analysis=true,
                        output_directory=analysis_directory,
                        save_everystep=true, dense=false,
                        l2=[0.05128549457329617], linf=[0.41441386613429104])

    analysis_path = joinpath(analysis_directory, "analysis.dat")
    history = accepted_step_history(analysis_path)
    @test SciMLBase.successful_retcode(sol)
    @test sol.stats.nreject > 0
    @test history.steps == collect(1:sol.stats.naccept)
    @test eltype(history.steps) == Int
    @test history.times == sol.t[2:end]
    @test history.dts≈diff(sol.t) rtol=1.0e-12 atol=eps(1.0)

    full_history = accepted_step_history(analysis_path; include_initial = true)
    @test full_history.steps == collect(0:sol.stats.naccept)
    @test full_history.times == sol.t
    @test first(full_history.dts) == 1.0e-2
    @test full_history.dts[2:end] == history.dts

    times, biases = mass_bias_history(analysis_path)
    @test times == sol.t
    # The uniform initial profile gives an independent reference for total storage.
    column_length = last(problem.domain)[1] - first(problem.domain)[1]
    initial_storage = water_content(-0.615, problem.equations) * column_length
    @test last(biases)≈mass_bias(last(sol.u), semi, initial_storage) atol=1.0e-14

    # Preserve the initial row and signed biases in the conservation-study schema.
    study_path = joinpath(analysis_directory, "conservation.dat")
    open(study_path, "w") do io
        println(io, "\n#accepted_step time_s dt_s mass_bias_m")
        println(io, "\n# Accepted-step samples")
        for row in zip(full_history.steps, full_history.times, full_history.dts, biases)
            println(io, join(row, ' '))
        end
    end
    columns = (; step_column = "accepted_step", time_column = "time_s", dt_column = "dt_s")
    @test accepted_step_history(study_path; columns...) == history
    @test accepted_step_history(study_path; columns..., include_initial = true) ==
          full_history
    @test mass_bias_history(study_path; time_column = "time_s",
                            mass_balance_column = "mass_bias_m") == (times, biases)
end

@trixi_testset "elixir_richards_celia_haverkamp.jl saved AMR meshes" begin
    elixir = joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_celia_haverkamp.jl")
    @test_trixi_include(elixir, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr=true, amr_interval=1, base_level=1, max_level=4,
                        save_mesh_history=true,
                        saveat=0.0:0.25:1.0, dense=false,
                        l2=[0.05128549457329617], linf=[0.41441386613429104])
    @test result.sol === sol
    @test sol.t == collect(0.0:0.25:1.0)
    @test length(result.mesh_history) == length(sol.u)
    @test length(first(result.mesh_history)) == 5
    @test length(last(result.mesh_history)) == 6
    @test result.mesh_history[2] === result.mesh_history[3]

    # The initial uniform profile is extracted on its original mesh after refinement.
    initial = solution_data_1d(result; index = 1, component = 2)
    initial_water = solution_data_1d(result; index = 1, component = 1)
    @test initial.time≈first(sol.t) atol=eps(1.0)
    @test initial.mesh_vertices_x≈first(result.mesh_history) rtol=1.0e-14 atol=eps(1.0)
    @test all(x -> first(problem.domain)[1] <= x <= last(problem.domain)[1],
              filter(isfinite, initial.x))
    @test filter(isfinite, initial.values) ≈ fill(-0.615, count(isfinite, initial.values))
    @test filter(isfinite, initial_water.values) ≈
          fill(water_content(-0.615, problem.equations),
               count(isfinite, initial_water.values))
    @test length(initial.x) == length(initial.values)
    @test all(edge -> count(==(edge), initial.x) == 2, initial.mesh_vertices_x[2:end-1])
    explicit = solution_data_1d(sol; index = 1, component = 2,
                                mesh_history = result.mesh_history)
    @test initial.time≈explicit.time rtol=1.0e-14 atol=eps(1.0)
    @test initial.x≈explicit.x rtol=1.0e-14 atol=eps(1.0)
    @test initial.values≈explicit.values rtol=1.0e-14 atol=eps(1.0)
    @test initial.mesh_vertices_x≈explicit.mesh_vertices_x rtol=1.0e-14 atol=eps(1.0)
    final = solution_data_1d(result; component = 2)
    @test final.time≈last(sol.t) rtol=1.0e-14 atol=eps(1.0)
    @test final.mesh_vertices_x≈last(result.mesh_history) rtol=1.0e-14 atol=eps(1.0)
    @test length(initial.x) != length(final.x)

    # Mesh-edge edits to a snapshot preserve shared saved meshes.
    shared_vertices = copy(result.mesh_history[2])
    snapshot = solution_data_1d(result; index = 2, component = 2)
    snapshot.mesh_vertices_x .*= 100
    @test result.mesh_history[2] == shared_vertices
    @test solution_data_1d(result; index = 3, component = 2).mesh_vertices_x ==
          shared_vertices

    @test_trixi_include(elixir, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr=true, amr_interval=1, base_level=1, max_level=4,
                        save_mesh_history=true,
                        l2=[0.05128549457329617], linf=[0.41441386613429104])
    @test length(sol.u) == 2
    @test length(result.mesh_history) == 2
    @test length(first(result.mesh_history)) == 5
    @test length(last(result.mesh_history)) == 6
    @test sol.t == [0.0, 1.0]

    # Solver defaults retain only the requested times for explicit save grids.
    for (requested_times, expected_times) in ((0.25, collect(0.0:0.25:1.0)),
                                              ([0.25, 0.75], [0.25, 0.75]),
                                              ([0.5, 1.0], [0.5, 1.0]))
        @test_trixi_include(elixir, tspan=(0.0, 1.0), initial_refinement_level=2,
                            amr=true, amr_interval=1, base_level=1, max_level=4,
                            save_mesh_history=true, saveat=requested_times, dense=false)
        @test SciMLBase.successful_retcode(sol)
        @test sol.t == expected_times
        @test length(result.mesh_history) == length(sol.u)
        for index in eachindex(sol.t)
            snapshot = solution_data_1d(result; index, component = 2)
            @test snapshot.time == sol.t[index]
            @test snapshot.mesh_vertices_x == result.mesh_history[index]
        end
    end

    @test_trixi_include(elixir, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr=true, amr_interval=1, base_level=1, max_level=4,
                        save_mesh_history=true,
                        save_everystep=true,
                        l2=[0.05128549457329617], linf=[0.41441386613429104])
    @test length(sol.u) == sol.stats.naccept + 1
    @test length(result.mesh_history) == length(sol.u)
    @test length(first(result.mesh_history)) == 5
    @test length(last(result.mesh_history)) == 6

    @test_trixi_include(elixir, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr=true, amr_interval=1, base_level=1, max_level=4,
                        save_mesh_history=true,
                        solve_options=(; save_start=false))
    @test length(sol.u) == 1
    @test length(result.mesh_history) == 1
    @test length(only(result.mesh_history)) == 6

    @test_trixi_include(elixir, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr=true, amr_interval=1, base_level=1, max_level=4,
                        save_mesh_history=true,
                        solve_options=(; save_end=false))
    @test length(sol.u) == 1
    @test length(result.mesh_history) == 1
    @test length(only(result.mesh_history)) == 5

    @test_trixi_include(elixir, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr=false, save_mesh_history=true,
                        saveat=0.0:0.25:1.0, dense=false)
    @test length(result.mesh_history) == length(sol.u)
    @test all(mesh -> mesh === first(result.mesh_history), result.mesh_history)
    # Static-mesh extraction agrees with and without a recorded mesh history.
    recorded = solution_data_1d(result; index = 1, component = 2)
    current = solution_data_1d(sol; index = 1, component = 2)
    # Reconstruction and live-mesh extraction may differ by floating-point roundoff.
    @test recorded.time≈current.time rtol=1.0e-14 atol=eps(1.0)
    @test recorded.x≈current.x rtol=1.0e-14 atol=eps(1.0)
    @test recorded.values≈current.values rtol=1.0e-14 atol=eps(1.0)
    @test recorded.mesh_vertices_x≈current.mesh_vertices_x rtol=1.0e-14 atol=eps(1.0)
end

@trixi_testset "elixir_richards_celia_haverkamp.jl normalized saturation indicator" begin
    # Normalization removes the affine change from water content to saturation.
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        amr=true, variable=effective_saturation,
                        tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        l2=[0.05128549457329617], linf=[0.41441386613429104])
end

@trixi_testset "elixir_richards_celia_new_mexico.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_new_mexico.jl"),
                        amr=false, initial_refinement_level=5,
                        l2=[6.561462104122631], linf=[9.250000021161606])

    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    @test extrema(Trixi.current_element_levels(mesh, dg, cache)) == (5, 5)
end

@trixi_testset "elixir_richards_celia_new_mexico.jl normalized saturation indicator" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_new_mexico.jl"),
                        amr=true, variable=effective_saturation,
                        tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        l2=[0.565135086335261], linf=[7.801851239219244])

    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    @test extrema(Trixi.current_element_levels(mesh, dg, cache)) == (1, 4)
end

@trixi_testset "elixir_richards_celia_haverkamp.jl pressure-head mapped error AMR" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        form=PressureHeadForm(),
                        error_control_block=state_variable_block,
                        error_control_mapping=water_content,
                        amr=true, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        l2=[0.05128898671388761], linf=[0.41441316902984093])

    @test SciMLBase.successful_retcode(sol)
    @test sol.stats.nreject > 0
    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    levels = Trixi.current_element_levels(mesh, dg, cache)
    @test extrema(levels) == (1, 4)
    @test length(sol.u[end]) > length(ode.u0)
end

@trixi_testset "elixir_richards_manufactured_solution.jl mixed form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        problem=HydrologicProblemRichardsManufacturedSolution(
                            boundary_conditions = :dirichlet_dirichlet),
                        l2=[5.842386475768436e-5], linf=[0.0003809050528057467])
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        problem=HydrologicProblemRichardsManufacturedSolution(
                            boundary_conditions = nothing),
                        l2=[5.842386475768436e-5], linf=[0.0003809050528057467])
end

@trixi_testset "elixir_richards_manufactured_solution.jl state-variable block" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        error_control_block=state_variable_block,
                        l2=[5.8423864976671844e-5], linf=[0.0003809050237838507])
end

@trixi_testset "elixir_richards_manufactured_solution.jl zero penalty factor" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        problem=HydrologicProblemRichardsManufacturedSolution(penalty_factor = 0),
                        l2=[9.188627055911749e-5], linf=[0.0005052944044751373])
end

@trixi_testset "elixir_richards_manufactured_solution.jl finite-diff Jacobian" begin
    import OrdinaryDiffEqRosenbrock

    # Retain support for graph-coloured finite-difference Jacobians
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        algorithm=default_algorithm(ode;
                                                    autodiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff()),
                        l2=[5.842386480880799e-5],
                        linf=[0.0003809050445998108])
end

@trixi_testset "elixir_richards_manufactured_solution.jl pressure-head form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        form=PressureHeadForm(), l2=[5.842386481713521e-5],
                        linf=[0.00038090499974052783])
end

@trixi_testset "elixir_richards_manufactured_solution.jl mapped error control" begin
    import OrdinaryDiffEqRosenbrock

    # Cover both final-stage estimators and weighted estimators such as Rodas5Pe.
    for (integration_algorithm, l2_reference, linf_reference) in
        ((OrdinaryDiffEqRosenbrock.Rodas4(), 0.006910145259419207, 0.025906733031605955),
         (OrdinaryDiffEqRosenbrock.Rodas42(), 0.006910141884804654, 0.025906634553169217),
         (OrdinaryDiffEqRosenbrock.Rodas4P(), 0.006910121384835785, 0.025906665817691632),
         (OrdinaryDiffEqRosenbrock.Rodas4P2(), 0.006910126223905007, 0.025906674745248104),
         (OrdinaryDiffEqRosenbrock.Rodas5(), 0.0069100575530365626, 0.025906469159871826),
         (OrdinaryDiffEqRosenbrock.Rodas5P(), 0.0069100589580171606, 0.025906458238520114),
         (OrdinaryDiffEqRosenbrock.Rodas5Pe(), 0.006910045004966427, 0.025906403901456265),
         (OrdinaryDiffEqRosenbrock.Rodas6P(), 0.006910058874553538, 0.02590647485709885))
        @testset "$(nameof(typeof(integration_algorithm)))" begin
            @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                         "elixir_richards_manufactured_solution.jl"),
                                algorithm=integration_algorithm, form=PressureHeadForm(),
                                error_control_block=state_variable_block,
                                error_control_mapping=water_content,
                                tspan=(0.0, 10.0), initial_refinement_level=2,
                                dt=0.5, reltol=1.0e-4, abstol=1.0e-8,
                                l2=[l2_reference], linf=[linf_reference])
            @test SciMLBase.successful_retcode(sol)
        end
    end
end

@trixi_testset "elixir_richards_manufactured_solution.jl mixed mapped error control" begin
    import OrdinaryDiffEqCore
    import OrdinaryDiffEqRosenbrock

    mapping_calls = Ref(0)
    function checked_water_content(value, equations)
        # Explicit mappings receive pressure head, including for the mixed form.
        @assert value < 0
        mapping_calls[] += 1
        return water_content(value, equations)
    end

    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        error_control_block=state_variable_block,
                        error_control_mapping=checked_water_content,
                        tspan=(0.0, 10.0), initial_refinement_level=2,
                        dt=0.5, reltol=1.0e-4, abstol=fill(1.0e-8, length(ode.u0)),
                        l2=[0.006910077076977245], linf=[0.025906520810020095])
    @test SciMLBase.successful_retcode(sol)
    @test mapping_calls[] > 0

    # Reinitialization must reset the wrapped PI controller's adaptive history.
    controller = HydroTrixi.StepsizeControllerMappedError(
        OrdinaryDiffEqCore.PIController(0.14, 0.08), state_variable_block,
        checked_water_content)
    integrator = SciMLBase.init(ode, default_algorithm(ode);
                                controller, dt = 0.5, reltol = 1.0e-4,
                                abstol = fill(1.0e-8, length(ode.u0)),
                                internalnorm = HydroTrixi.variable_block_norm(
                                    evolved_variable_block, semi),
                                save_everystep = false)
    SciMLBase.step!(integrator)
    SciMLBase.reinit!(integrator)
    reinitialized_solution = SciMLBase.solve!(integrator)
    @test SciMLBase.successful_retcode(reinitialized_solution)
    @test last(reinitialized_solution.u)≈last(sol.u) rtol=1.0e-10

    options = (; dt = 0.5, error_control_block = state_variable_block,
               error_control_mapping = checked_water_content)
    for unsupported in (OrdinaryDiffEqRosenbrock.Rosenbrock23(),
                        OrdinaryDiffEqRosenbrock.Rodas5Pr())
        @test_throws ArgumentError solve_implicit(ode, unsupported; options...)
    end
    @test_throws ArgumentError solve_implicit(ode; options...,
                                              reltol = fill(1.0e-4, length(ode.u0)))
    @test_throws ArgumentError solve_implicit(ode; options...,
                                              step_limiter = (u, integrator, p, t) -> nothing)
    # Fixed stepping does not evaluate the optional mapping.
    mapping_calls[] = 0
    fixed = solve_implicit(ode; options..., adaptive = false).sol
    @test SciMLBase.successful_retcode(fixed)
    @test mapping_calls[] == 0
end

@testset "Richards manufactured solution Dirichlet-Neumann convergence" begin
    preset = HydrologicProblemRichardsManufacturedSolution(
        boundary_conditions = :dirichlet_neumann)
    # Retain the custom unpenalized Dirichlet condition used by this regression.
    boundaries = (; x_neg = Trixi.BoundaryConditionDirichlet(preset.initial_condition),
                   x_pos = preset.boundary_conditions.x_pos)
    custom = HydrologicProblemRichardsManufacturedSolution(boundary_conditions = boundaries)

    for problem in (custom, preset)
        eocs, _ = Trixi.convergence_test(@__MODULE__,
                                         joinpath(EXAMPLES_DIR, "elixirs",
                                                  "elixir_richards_manufactured_solution.jl"),
                                         3; problem = problem,
                                         initial_refinement_level = 4)
        mean_convergence = Trixi.calc_mean_convergence(eocs)
        @test isapprox(mean_convergence[:l2], [4.0]; rtol = 0.1)
        @test isapprox(mean_convergence[:linf], [4.0]; rtol = 0.1)
    end
end

@testset "elixir_richards_celia_haverkamp.jl AMR mass bias" begin
    analysis_directory = mktempdir()
    Trixi.trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl");
                        tspan = (0.0, 1.0), amr = true, run_simulation = false,
                        save_analysis = true, output_directory = analysis_directory)

    # AMR must not trigger DAE reinitialization when it leaves the mesh unchanged.
    sol = solve(ode, default_algorithm(ode); dt = 1.0e-2, adaptive = true,
                initializealg = SciMLBase.CheckInit(),
                reltol = 1.0e-7, abstol = 1.0e-11, save_everystep = false,
                save_start = false,
                save_end = true, maxiters = typemax(Int), callback = callbacks)
    analysis_path = joinpath(analysis_directory, "analysis.dat")
    _, biases = HydroTrixi.mass_bias_history(analysis_path)
    @test abs(last(biases)) < 1.0e-12
end

@testset "elixir_richards_celia_haverkamp.jl AMR Jacobian" begin
    # Scheduled AMR callback to keep simulation topologies consistent between runs
    function solve_scheduled_amr(ode, semi, mesh, amr_callback, adaptation_times;
                                 algorithm = default_algorithm(ode))
        topology_history = Tuple{Float64, Vector{Int}}[]
        scheduled_times = Set(adaptation_times)
        condition = (u, t, integrator) -> t in scheduled_times
        affect! = function (integrator)
            amr_callback.affect!(integrator)
            # Both refinement and coarsening keep the reused spatial storage consistent.
            cache_parabolic = semi.cache_parabolic
            n_elements = Trixi.nelements(semi.semi_base.solver, semi.semi_base.cache)
            @test size(cache_parabolic.parabolic_container.u_transformed, 3) == n_elements
            @test size(cache_parabolic.parabolic_boundaries.flux_values) ==
                  (2, Trixi.nvariables(semi.semi_base),
                   Trixi.nboundaries(semi.semi_base.cache.boundaries))
            push!(topology_history,
                  (integrator.t, copy(Trixi.leaf_cells(mesh.tree))))
            return nothing
        end
        callback = DiscreteCallback(condition, affect!;
                                    save_positions = (false, false))

        # Use the same time grid for dense and sparse Jacobian runs
        solution = solve(ode, algorithm;
                         dt = 1.0e-2, adaptive = false, saveat = Float64[],
                         Trixi.ode_default_options()...,
                         callback = callback,
                         tstops = adaptation_times)
        return (; solution, topology_history)
    end

    elixir = joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_celia_haverkamp.jl")
    adaptation_times = collect(30.0:30.0:330.0)

    # Test dense and sparse Jacobian runs for both mixed and pressure-head forms
    for form in (MixedForm(), PressureHeadForm())
        @testset "$(nameof(typeof(form)))" begin
            dense, sparse = map((DenseJacobian(), SparseJacobian())) do jacobian_strategy
                # Bound the dense Jacobian size while exercising normalized AMR.
                Trixi.trixi_include(@__MODULE__, elixir;
                                    form = form, jacobian = jacobian_strategy, amr = true,
                                    initial_refinement_level = 2, max_level = 6,
                                    run_simulation = false)
                solve_scheduled_amr(ode, semi, mesh, amr_callback, adaptation_times)
            end

            @test successful_retcode(dense.solution)
            @test successful_retcode(sparse.solution)
            @test length(sparse.topology_history) == length(adaptation_times)
            @test dense.topology_history == sparse.topology_history
            @test maximum(abs,
                          last(dense.solution.u) .- last(sparse.solution.u)) < 1.0e-9
        end
    end
end

@testset "elixir_richards_closed_column.jl mass conservation" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_closed_column.jl"),
                        saveat=0.0:10.0:360.0, amr=true,
                        l2=[0.0002775234378470748], linf=[0.0005159084431309857])

    storage = [HydroTrixi.water_content_integral(u_ode, semi) for u_ode in sol.u]
    initial_storage = first(storage)
    @test isapprox(storage, fill(initial_storage, length(storage));
                   rtol = 0, atol = 100 * eps(initial_storage))
end
