using Test
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
