using Test
using HydroTrixi
using Trixi: trixi_include
import Trixi
using TrixiTest

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")

macro test_trixi_include(args...)
    esc(Expr(:macrocall, Symbol("@test_trixi_include_base"), __source__, args...))
end

@trixi_testset "elixir_diffusion_1d_dirichlet_dirichlet.jl" begin
    @test_trixi_include joinpath(EXAMPLES_DIR,
                                 "elixir_diffusion_1d_dirichlet_dirichlet.jl")
end

@trixi_testset "elixir_diffusion_1d_mixed_dirichlet_neumann.jl" begin
    @test_trixi_include joinpath(EXAMPLES_DIR,
                                 "elixir_diffusion_1d_mixed_dirichlet_neumann.jl")
end

@trixi_testset "elixir_richards_celia_1990.jl" begin
    @test_trixi_include joinpath(EXAMPLES_DIR,
                                 "elixir_richards_celia_1990.jl")
end

@testset "elixir_richards_closed_column.jl mass conservation" begin
    trixi_include(joinpath(EXAMPLES_DIR, "elixir_richards_closed_column.jl");
                  final_time = 360.0,
                  saveat = 0.0:10.0:360.0)

    storage = similar(sol.t, length(sol.t))
    rate = similar(sol.t, length(sol.t))
    for (i, (u_ode, t)) in enumerate(zip(sol.u, sol.t))
        du_ode = similar(u_ode)
        Trixi.default_rhs(semi)(du_ode, u_ode, semi, t)
        storage[i] = Trixi.analyze(HydroTrixi.water_content, du_ode, u_ode, t, semi)
        rate[i] = Trixi.analyze(HydroTrixi.water_content_timederivative, du_ode, u_ode, t,
                                semi)
    end

    initial_storage = first(storage)
    relative_storage_drift = maximum(abs.(storage .- initial_storage)) /
                             abs(initial_storage)
    relative_rate = maximum(abs.(rate)) / abs(initial_storage)

    @test relative_storage_drift < 5.0e-6
    @test relative_rate < 5.0e-6

    base_semi = semi.semi_base
    u_ode = sol.u[end]
    t = sol.t[end]

    storage_direct = Trixi.integrate(Trixi.cons2cons,
                                     HydroTrixi.evolved_variable_view(u_ode),
                                     base_semi;
                                     normalize = false)
    @test HydroTrixi.evolved_variables_integral(u_ode, semi) ≈ storage_direct
    @test Trixi.analyze(HydroTrixi.water_content, u_ode, u_ode, t, semi) ≈
          first(storage_direct)

    du_ode = similar(u_ode)
    Trixi.default_rhs(semi)(du_ode, u_ode, semi, t)

    rate_direct = Trixi.integrate(Trixi.cons2cons,
                                  HydroTrixi.evolved_variable_view(du_ode),
                                  base_semi;
                                  normalize = false)
    @test HydroTrixi.evolved_variables_timederivative(u_ode, semi, t) ≈ rate_direct
    @test Trixi.analyze(HydroTrixi.water_content_timederivative, du_ode, u_ode, t, semi) ≈
          first(rate_direct)

    @test Trixi.analyze(HydroTrixi.water_content, du_ode, u_ode, t, semi) ≈
          first(storage_direct)
end

if isempty(ARGS) || "visualization" in ARGS
    include("visualization_extension.jl")
end
