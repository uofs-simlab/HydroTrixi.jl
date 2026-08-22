using Test
using HydroTrixi
import LinearSolve
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
    end
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl" begin
    import LinearSolve
    import OrdinaryDiffEqRosenbrock
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        l2=[4.688250908054879e-5], linf=[0.00035212174570349586])
    @test sol.alg.autodiff isa OrdinaryDiffEqRosenbrock.AutoForwardDiff
    @test sol.alg.linsolve isa LinearSolve.LUFactorization
    @test sol.alg.max_jac_age == 1
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl implicit" begin
    import LinearSolve
    import OrdinaryDiffEqRosenbrock
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        algorithm=default_algorithm(ode), jacobian=SparseJacobian(),
                        dt=1.0e-2, adaptive=true,
                        reltol=1.0e-9, abstol=1.0e-11,
                        l2=[4.688250908054879e-5], linf=[0.00035212174570349586])
    @test sol.alg.autodiff isa OrdinaryDiffEqRosenbrock.AutoSparse
    @test sol.alg.autodiff.dense_ad isa OrdinaryDiffEqRosenbrock.AutoForwardDiff
    @test sol.alg.linsolve isa LinearSolve.KLUFactorization
    @test sol.alg.max_jac_age == 1
end

@trixi_testset "elixir_diffusion_1d_mixed_dirichlet_neumann.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_mixed_dirichlet_neumann.jl"),
                        l2=[2.7083226488116088e-5], linf=[0.00022679747793086236])
end

@trixi_testset "elixir_richards_celia_1990.jl" begin
    @test_trixi_include joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_celia_1990.jl")
end

@trixi_testset "elixir_richards_manufactured_solution.jl mixed form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        l2=[4.0696092211162146e-5], linf=[0.0003809050528035818])
end

@testset "elixir_richards_manufactured_solution.jl zero penalty factor" begin
    problem = HydrologicProblemRichardsManufacturedSolution(penalty_factor = 0)
    Trixi.trixi_include(@__MODULE__,
                        joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl");
                        problem = problem)
    errors = analysis_callback(sol)
    @test successful_retcode(sol)
    @test errors.l2≈[6.174720607183763e-5] rtol=1.0e-10
    @test errors.linf≈[0.0005052944044764973] rtol=1.0e-10
end

@testset "elixir_richards_manufactured_solution.jl finite-diff Jacobian" begin
    # Retain support for graph-coloured finite-difference Jacobians
    autodiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff()
    finite_diff_algorithm = default_algorithm(ode; autodiff)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        algorithm=finite_diff_algorithm, l2=[4.0696092224417466e-5],
                        linf=[0.0003809050439483319])
end

@trixi_testset "elixir_richards_manufactured_solution.jl pressure-head form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        form=PressureHeadForm(), l2=[4.069609136528145e-5],
                        linf=[0.000380904999743803])
end

@testset "state-variable norm" begin
    elixir = joinpath(EXAMPLES_DIR, "elixirs",
                      "elixir_richards_manufactured_solution.jl")

    for form in (PressureHeadForm(), MixedForm())
        Trixi.trixi_include(@__MODULE__, elixir; form = form, run_simulation = false)
        error_norm = state_variable_norm(semi)
        n_state_variables = length(HydroTrixi.state_variable_view(ode.u0, semi))

        test_state = if form isa PressureHeadForm
            fill(-2.0, n_state_variables)
        else
            vcat(fill(100.0, n_state_variables), fill(-2.0, n_state_variables))
        end

        @test @inferred(error_norm(-2.0, 0.0)) == 2.0
        @test @inferred(error_norm(test_state, 0.0)) == 2.0
    end

    semi_standard = SemidiscretizationImplicit(semi.semi_base,
                                               TemporalOperatorStandard())
    standard_norm = state_variable_norm(semi_standard)
    @test @inferred(standard_norm(fill(-2.0, length(ode.u0) ÷ 2), 0.0)) == 2.0

    celia_elixir = joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_celia_1990.jl")
    Trixi.trixi_include(@__MODULE__, celia_elixir; run_simulation = false)
    error_norm = state_variable_norm(semi)
    n_state_variables = length(HydroTrixi.state_variable_view(ode.u0, semi))
    test_state = vcat(fill(100.0, n_state_variables), fill(-2.0, n_state_variables),
                      [1.0e6, -1.0e6])
    @test @inferred(error_norm(test_state, 0.0)) == 2.0
end

@testset "evolved-variable norm" begin
    elixir = joinpath(EXAMPLES_DIR, "elixirs",
                      "elixir_richards_manufactured_solution.jl")

    for form in (PressureHeadForm(), MixedForm())
        Trixi.trixi_include(@__MODULE__, elixir; form = form, run_simulation = false)
        error_norm = evolved_variable_norm(semi)
        n_evolved_variables = length(HydroTrixi.evolved_variable_view(ode.u0, semi))

        test_state = if form isa PressureHeadForm
            fill(-2.0, n_evolved_variables)
        else
            vcat(fill(-2.0, n_evolved_variables), fill(100.0, n_evolved_variables))
        end

        @test @inferred(error_norm(-2.0, 0.0)) == 2.0
        @test @inferred(error_norm(test_state, 0.0)) == 2.0
    end

    semi_standard = SemidiscretizationImplicit(semi.semi_base,
                                               TemporalOperatorStandard())
    standard_norm = evolved_variable_norm(semi_standard)
    @test @inferred(standard_norm(fill(-2.0, length(ode.u0) ÷ 2), 0.0)) == 2.0

    celia_elixir = joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_celia_1990.jl")
    Trixi.trixi_include(@__MODULE__, celia_elixir; run_simulation = false)
    error_norm = evolved_variable_norm(semi)
    n_evolved_variables = length(HydroTrixi.evolved_variable_view(ode.u0, semi))
    test_state = vcat(fill(-2.0, n_evolved_variables),
                      fill(100.0, n_evolved_variables), [1.0e6, -1.0e6])
    @test @inferred(error_norm(test_state, 0.0)) == 2.0
end

@testset "elixir_richards_celia_1990.jl AMR mass bias" begin
    Trixi.trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_1990.jl");
                        tspan = (0.0, 1.0), amr = true, run_simulation = false)

    sol = solve(ode, default_algorithm(ode); dt = 1.0e-2, adaptive = true,
                reltol = 1.0e-7, abstol = 1.0e-11, save_everystep = false,
                save_start = false,
                save_end = true, maxiters = typemax(Int), callback = callbacks)
    @test abs(HydroTrixi.mass_bias(sol.u[end], semi)) < 1.0e-12
end

@testset "elixir_richards_celia_1990.jl AMR Jacobian" begin
    # Scheduled AMR callback to keep simulation topologies consistent between runs
    function solve_scheduled_amr(ode, semi, mesh, amr_callback, adaptation_times;
                                 algorithm = default_algorithm(ode))
        topology_history = Tuple{Float64, Vector{Int}}[]
        scheduled_times = Set(adaptation_times)
        condition = (u, t, integrator) -> t in scheduled_times
        affect! = function (integrator)
            amr_callback.affect!(integrator)
            push!(topology_history,
                  (integrator.t, copy(Trixi.leaf_cells(mesh.tree))))
            return nothing
        end
        callback = DiscreteCallback(condition, affect!;
                                    save_positions = (false, false))

        solution = solve(ode, algorithm;
                         dt = 1.0e-2, adaptive = true, reltol = 1.0e-7,
                         abstol = 1.0e-11, saveat = Float64[],
                         Trixi.ode_default_options()...,
                         callback = callback,
                         tstops = adaptation_times)
        return (; solution, topology_history)
    end

    elixir = joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_celia_1990.jl")
    adaptation_times = collect(30.0:30.0:330.0)

    # Test dense and sparse Jacobian runs for both mixed and pressure-head forms
    for form in (MixedForm(), PressureHeadForm())
        @testset "$(nameof(typeof(form)))" begin
            dense, sparse = map((DenseJacobian(), SparseJacobian())) do jacobian_strategy
                Trixi.trixi_include(@__MODULE__, elixir;
                                    form = form, jacobian = jacobian_strategy, amr = true,
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
    trixi_include(joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_closed_column.jl");
                  saveat = 0.0:10.0:360.0)

    storage = [only(HydroTrixi.evolved_variables_integral(u_ode, semi)) for u_ode in sol.u]
    initial_storage = first(storage)
    @test isapprox(storage, fill(initial_storage, length(storage));
                   rtol = 0, atol = 100 * eps(initial_storage))
end

if isempty(ARGS) || "visualization" in ARGS
    include("visualization_extension.jl")
end
