using CairoMakie
using LaTeXStrings

@testset "visualization extension smoke test" begin
    trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                           "elixir_diffusion_1d_mixed_dirichlet_neumann.jl"))

    mktempdir() do tmpdir
        output_path = joinpath(tmpdir, "diffusion_solution.png")
        fig = plot_solution_1d(sol; exact_solution = exact_solution,
                               output_path = output_path)

        @test fig isa CairoMakie.Figure
        @test isfile(output_path)
    end
end
