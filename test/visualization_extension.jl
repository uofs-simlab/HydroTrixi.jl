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

@testset "arbitrary convergence series" begin
    x = [8.0, 16.0, 32.0]
    series = ((; x, errors = (x .^ -2, 2 .* x .^ -2, x .^ -3),
               labels = (L"$L^2$", L"$L^\infty$", L"$H^1$"), color = 1,
               marker = :rect),)

    mktempdir() do tmpdir
        output_path = joinpath(tmpdir, "convergence.pdf")
        fig = plot_convergence_1d(series; output_path = output_path,
                                  triangle_order = 3)

        @test fig isa CairoMakie.Figure
        @test isfile(output_path)
    end
end
