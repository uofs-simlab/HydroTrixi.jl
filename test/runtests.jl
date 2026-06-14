using Test
using HydroTrixi
using SciMLBase: solve
using Trixi: trixi_include
import Trixi
using TrixiTest

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")

macro test_trixi_include(args...)
    esc(Expr(:macrocall, Symbol("@test_trixi_include_base"), __source__, args...))
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl"),
                        l2=[4.688250908054879e-5],
                        linf=[0.00035212174570349586])
end

@trixi_testset "elixir_diffusion_1d_mixed_dirichlet_neumann.jl" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_diffusion_1d_mixed_dirichlet_neumann.jl"),
                        l2=[2.7083226488116088e-5],
                        linf=[0.00022679747793086236])
end

@trixi_testset "elixir_richards_celia_1990.jl" begin
    @test_trixi_include joinpath(EXAMPLES_DIR,
                                 "elixir_richards_celia_1990.jl")
end

@trixi_testset "elixir_richards_manufactured_solution.jl mixed form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_richards_manufactured_solution.jl"),
                        form=MixedForm(),
                        final_time=120.0,
                        initial_refinement_level=4,
                        polydeg=3,
                        reltol=1.0e-9,
                        abstol=1.0e-11,
                        saveat=Float64[],
                        l2=[6.174720620257904e-5],
                        linf=[0.0005052944044678376])
end

@trixi_testset "elixir_richards_manufactured_solution.jl pressure-head form" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "elixir_richards_manufactured_solution.jl"),
                        form=PressureHeadForm(),
                        final_time=120.0,
                        initial_refinement_level=4,
                        polydeg=3,
                        reltol=1.0e-9,
                        abstol=1.0e-11,
                        saveat=Float64[],
                        l2=[6.174720554384852e-5],
                        linf=[0.0005052944065825626])
end

@testset "elixir_richards_celia_1990_amr.jl mass bias" begin
    Trixi.trixi_include(joinpath(EXAMPLES_DIR, "elixir_richards_celia_1990_amr.jl");
                        final_time = 1.0,
                        run_simulation = false,
                        saveat = Float64[],
                        amr_interval = 30,
                        max_level = 6,
                        reltol = 1.0e-7,
                        abstol = 1.0e-11)

    sol = solve(ode, default_algorithm(semi);
                dt = dt,
                adaptive = adaptive,
                reltol = reltol,
                abstol = abstol,
                save_everystep = false,
                save_start = false,
                save_end = true,
                maxiters = typemax(Int),
                callback = callbacks)
    @test abs(HydroTrixi.mass_bias(sol.u[end], semi)) < 1.0e-12
end

@testset "elixir_richards_celia_1990_amr.jl pressure-head form" begin
    Trixi.trixi_include(joinpath(EXAMPLES_DIR, "elixir_richards_celia_1990_amr.jl");
                        final_time = 1.0,
                        form = PressureHeadForm(),
                        run_simulation = false,
                        saveat = Float64[],
                        amr_interval = 30,
                        max_level = 6,
                        reltol = 1.0e-7,
                        abstol = 1.0e-11)

    @test semi.operator_temporal isa HydroTrixi.TemporalOperatorCapacity

    pressure_head_solution = solve(ode, default_algorithm(semi);
                                   dt = dt,
                                   adaptive = adaptive,
                                   reltol = reltol,
                                   abstol = abstol,
                                   save_everystep = false,
                                   save_start = false,
                                   save_end = true,
                                   maxiters = typemax(Int),
                                   callback = callbacks)
    @test pressure_head_solution.t[end] ≈ 1.0
end

@testset "elixir_richards_closed_column.jl mass conservation" begin
    trixi_include(joinpath(EXAMPLES_DIR, "elixir_richards_closed_column.jl");
                  final_time = 360.0,
                  saveat = 0.0:10.0:360.0)

    storage = [only(HydroTrixi.evolved_variables_integral(u_ode, semi))
               for u_ode in sol.u]
    initial_storage = first(storage)
    @test isapprox(storage, fill(initial_storage, length(storage));
                   rtol = 0, atol = 100 * eps(initial_storage))
end

if isempty(ARGS) || "visualization" in ARGS
    include("visualization_extension.jl")
end
