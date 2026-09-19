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
    end
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl dense Jacobian" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        l2=[4.688250908054879e-5], linf=[0.00035212174570349586])
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl sparse Jacobian" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        algorithm=default_algorithm(ode), jacobian=SparseJacobian(),
                        dt=1.0e-2, adaptive=true,
                        reltol=1.0e-9, abstol=1.0e-11,
                        l2=[4.688250908054879e-5], linf=[0.00035212174570349586])
end

@trixi_testset "elixir_diffusion_1d_mixed_dirichlet_neumann.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_diffusion_1d_mixed_dirichlet_neumann.jl"),
                        l2=[2.7083226488116088e-5], linf=[0.00022679747793086236])
end

@trixi_testset "elixir_richards_celia_haverkamp.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        amr=true,
                        l2=[0.2312163774363683], linf=[0.4080000141652919])

    # Normalization remains finite for constant and nearly constant states.
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    u = fill(-0.5, size(cache.elements.node_coordinates))
    @test all(iszero, amr_indicator(u, mesh, equations, dg, cache))
    u[1] += eps()
    @test all(isfinite, amr_indicator(u, mesh, equations, dg, cache))
end

@trixi_testset "elixir_richards_celia_haverkamp.jl normalized saturation indicator" begin
    # Normalization removes the affine change from water content to saturation.
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        amr=true, variable=effective_saturation,
                        l2=[0.2312163774363683], linf=[0.4080000141652919])
end

@trixi_testset "elixir_richards_celia_new_mexico.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_new_mexico.jl"),
                        l2=[6.748850474885821], linf=[9.250000021161606])

    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    @test maximum(Trixi.current_element_levels(mesh, dg, cache)) == 10
end

@trixi_testset "elixir_richards_celia_new_mexico.jl normalized saturation indicator" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_new_mexico.jl"),
                        variable=effective_saturation,
                        l2=[6.748850474885821], linf=[9.250000021161606])

    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    @test maximum(Trixi.current_element_levels(mesh, dg, cache)) == 10
end

@trixi_testset "elixir_richards_celia_haverkamp.jl pressure-head mapped error AMR" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        form=PressureHeadForm(), error_control_variables=water_content,
                        amr=true, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, amr_base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        l2=[0.04980091115465496], linf=[0.41441316902984104])

    @test SciMLBase.successful_retcode(sol)
    @test sol.stats.nreject > 0
    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    levels = Trixi.current_element_levels(mesh, dg, cache)
    @test extrema(levels) == (1, 4)
    @test length(sol.u[end]) > length(ode.u0)
end

@trixi_testset "elixir_richards_celia_haverkamp.jl mixed evolved error AMR" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        error_control_variables=water_content,
                        amr=true, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, amr_base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        l2=[0.04979785267197109], linf=[0.41441386613429065])

    @test SciMLBase.successful_retcode(sol)
    @test sol.stats.naccept == 83
    @test sol.stats.nreject > 0
    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    levels = Trixi.current_element_levels(mesh, dg, cache)
    @test extrema(levels) == (1, 4)
    @test length(sol.u[end]) > length(ode.u0)
end

@trixi_testset "elixir_richards_manufactured_solution.jl mixed form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        l2=[4.0696092211162146e-5], linf=[0.0003809050528035818])
end

@trixi_testset "elixir_richards_manufactured_solution.jl state-variable norm" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        internalnorm=state_variable_norm(semi),
                        l2=[4.069609236532572e-5], linf=[0.00038090502378751445])
end

@trixi_testset "elixir_richards_manufactured_solution.jl zero penalty factor" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        problem=HydrologicProblemRichardsManufacturedSolution(penalty_factor = 0),
                        l2=[6.174720607183763e-5], linf=[0.0005052944044764973])
end

@trixi_testset "elixir_richards_manufactured_solution.jl finite-diff Jacobian" begin
    import OrdinaryDiffEqRosenbrock

    # Retain support for graph-coloured finite-difference Jacobians
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        algorithm=default_algorithm(ode;
                                                    autodiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff()),
                        l2=[4.0696092224417466e-5],
                        linf=[0.0003809050439483319])
end

@trixi_testset "elixir_richards_manufactured_solution.jl pressure-head form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_manufactured_solution.jl"),
                        form=PressureHeadForm(), l2=[4.069609136528145e-5],
                        linf=[0.000380904999743803])
end

@trixi_testset "elixir_richards_manufactured_solution.jl mapped error control" begin
    import OrdinaryDiffEqRosenbrock

    # Cover both final-stage estimators and weighted estimators such as Rodas5Pe.
    for (integration_algorithm, l2_reference, linf_reference) in
        ((OrdinaryDiffEqRosenbrock.Rodas4(), 0.004825104682560886, 0.025906733031606066),
         (OrdinaryDiffEqRosenbrock.Rodas42(), 0.004825116227347337, 0.025906634553169328),
         (OrdinaryDiffEqRosenbrock.Rodas4P(), 0.004825093476699578, 0.025906665817691854),
         (OrdinaryDiffEqRosenbrock.Rodas4P2(), 0.004825096265966346, 0.02590667474524755),
         (OrdinaryDiffEqRosenbrock.Rodas5(), 0.0048250484071016, 0.025906469159872048),
         (OrdinaryDiffEqRosenbrock.Rodas5P(), 0.004825051926237008, 0.025906458238520003),
         (OrdinaryDiffEqRosenbrock.Rodas5Pe(), 0.004825045396870383, 0.025906403901456265),
         (OrdinaryDiffEqRosenbrock.Rodas6P(), 0.004825050923678062, 0.025906474857098738))
        @testset "$(nameof(typeof(integration_algorithm)))" begin
            @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                         "elixir_richards_manufactured_solution.jl"),
                                algorithm=integration_algorithm, form=PressureHeadForm(),
                                error_control_variables=water_content,
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
                        error_control_variables=checked_water_content,
                        tspan=(0.0, 10.0), initial_refinement_level=2,
                        dt=0.5, reltol=1.0e-4, abstol=fill(1.0e-8, length(ode.u0)),
                        l2=[0.004825064501334634], linf=[0.025906520810019762])
    @test SciMLBase.successful_retcode(sol)
    @test mapping_calls[] > 0

    # Reinitialization must reset the wrapped PI controller's adaptive history.
    controller = HydroTrixi.StepsizeControllerMappedError(
        OrdinaryDiffEqCore.PIController(0.14, 0.08), checked_water_content)
    integrator = SciMLBase.init(ode, default_algorithm(ode);
                                controller, dt = 0.5, reltol = 1.0e-4,
                                abstol = fill(1.0e-8, length(ode.u0)),
                                internalnorm = evolved_variable_norm(semi),
                                save_everystep = false)
    SciMLBase.step!(integrator)
    SciMLBase.reinit!(integrator)
    reinitialized_solution = SciMLBase.solve!(integrator)
    @test SciMLBase.successful_retcode(reinitialized_solution)
    @test last(reinitialized_solution.u)≈last(sol.u) rtol=1.0e-10

    options = (; dt = 0.5, error_control_variables = checked_water_content)
    for unsupported in (OrdinaryDiffEqRosenbrock.Rosenbrock23(),
                        OrdinaryDiffEqRosenbrock.Rodas5Pr())
        @test_throws ArgumentError solve_implicit(ode, unsupported; options...)
    end
    @test_throws ArgumentError solve_implicit(ode; options...,
                                              reltol = fill(1.0e-4, length(ode.u0)))
    @test_throws ArgumentError solve_implicit(ode; options...,
                                              step_limiter = (u, integrator, p, t) -> nothing)

    # Fixed stepping does not evaluate the optional conversion.
    mapping_calls[] = 0
    fixed = solve_implicit(ode; options..., adaptive = false)
    @test SciMLBase.successful_retcode(fixed)
    @test mapping_calls[] == 0
end

@testset "Richards manufactured solution Dirichlet-Neumann convergence" begin
    base_problem = HydrologicProblemRichardsManufacturedSolution()
    dirichlet_left = Trixi.BoundaryConditionDirichlet(base_problem.initial_condition)
    right_flux = HydroTrixi.richards_manufactured_right_boundary_flux
    neumann_right = Trixi.BoundaryConditionNeumann(right_flux)
    boundary_conditions = (; x_neg = dirichlet_left, x_pos = neumann_right)
    problem = HydrologicProblem(equations = base_problem.equations,
                                state_to_evolved = base_problem.state_to_evolved,
                                evolved_to_state = base_problem.evolved_to_state,
                                initial_condition = base_problem.initial_condition,
                                boundary_conditions = boundary_conditions,
                                source_terms = base_problem.source_terms,
                                domain = base_problem.domain, tspan = base_problem.tspan)

    eocs, _ = Trixi.convergence_test(@__MODULE__,
                                     joinpath(EXAMPLES_DIR, "elixirs",
                                              "elixir_richards_manufactured_solution.jl"),
                                     3; problem = problem,
                                     initial_refinement_level = 4)
    mean_convergence = Trixi.calc_mean_convergence(eocs)
    @test isapprox(mean_convergence[:l2], [4.0]; rtol = 0.1)
    @test isapprox(mean_convergence[:linf], [4.0]; rtol = 0.1)
end

@testset "elixir_richards_celia_haverkamp.jl AMR mass bias" begin
    Trixi.trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl");
                        tspan = (0.0, 1.0), amr = true, run_simulation = false)

    # AMR must not trigger DAE reinitialization when it leaves the mesh unchanged.
    sol = solve(ode, default_algorithm(ode); dt = 1.0e-2, adaptive = true,
                initializealg = SciMLBase.CheckInit(),
                reltol = 1.0e-7, abstol = 1.0e-11, save_everystep = false,
                save_start = false,
                save_end = true, maxiters = typemax(Int), callback = callbacks)
    @test abs(HydroTrixi.mass_bias(sol.u[end], semi)) < 1.0e-12
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

    elixir = joinpath(EXAMPLES_DIR, "elixirs", "elixir_richards_celia_haverkamp.jl")
    adaptation_times = collect(30.0:30.0:330.0)

    # Test dense and sparse Jacobian runs for both mixed and pressure-head forms
    for form in (MixedForm(), PressureHeadForm())
        @testset "$(nameof(typeof(form)))" begin
            dense, sparse = map((DenseJacobian(), SparseJacobian())) do jacobian_strategy
                # Bound the dense Jacobian size while exercising normalized AMR.
                Trixi.trixi_include(@__MODULE__, elixir;
                                    form = form, jacobian = jacobian_strategy, amr = true,
                                    max_level = 6,
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
                        l2=[0.00027752343492303905], linf=[0.000515908443132318])

    storage = [only(HydroTrixi.evolved_variables_integral(u_ode, semi)) for u_ode in sol.u]
    initial_storage = first(storage)
    @test isapprox(storage, fill(initial_storage, length(storage));
                   rtol = 0, atol = 100 * eps(initial_storage))
end

if isempty(ARGS) || "visualization" in ARGS
    include("visualization_extension.jl")
end
