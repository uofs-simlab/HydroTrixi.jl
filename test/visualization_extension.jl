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
        data = solution_data_1d(sol)
        points = fig.content[1].scene.plots[1][1][]
        @test first.(points)≈data.x nans=true rtol=1.0e-6
        @test last.(points)≈data.values nans=true rtol=1.0e-6

        # Extracted profiles use the same plotting path as ODE solutions.
        array_path = joinpath(tmpdir, "diffusion_solution_arrays.png")
        array_fig = plot_solution_1d(data; exact_solution, output_path = array_path)
        @test array_fig isa CairoMakie.Figure
        @test isfile(array_path)
        array_points = array_fig.content[1].scene.plots[1][1][]
        @test first.(array_points)≈first.(points) nans=true rtol=1.0e-6
        @test last.(array_points)≈last.(points) nans=true rtol=1.0e-6

        # The same public method accepts independently supplied arrays and mesh guides.
        profile = (; time = 0.5, x = [0.0, 0.1, 0.2], values = [-0.6, -0.5, -0.4],
                    mesh_vertices_x = [0.0, 0.2])
        profile_path = joinpath(tmpdir, "profile.pdf")
        profile_fig = plot_solution_1d(profile; output_path = profile_path,
                                        xlabel = "Depth (m)", ylabel = "Pressure head (m)",
                                        ylims = (-0.7, -0.3), show_element_boundaries = true)
        @test profile_fig isa CairoMakie.Figure
        @test isfile(profile_path)
        profile_axis = profile_fig.content[1]
        @test profile_axis.xlabel[] == "Depth (m)"
        @test profile_axis.ylabel[] == "Pressure head (m)"
        @test profile_axis.limits[][2] == (-0.7, -0.3)
        @test length(profile_axis.scene.plots) == 2
        @test profile_axis.limits[][1] == (0.0, 0.2)
    end
end

@trixi_testset "elixir_richards_celia_haverkamp.jl time-step plots" begin
    import CairoMakie

    analysis_directory = mktempdir()
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixirs",
                                 "elixir_richards_celia_haverkamp.jl"),
                        amr=true, tspan=(0.0, 1.0), initial_refinement_level=2,
                        amr_interval=1, base_level=1, max_level=4,
                        adapt_initial_condition=false,
                        analysis_interval=1, save_analysis=true,
                        output_directory=analysis_directory,
                        l2=[0.05128549457329617], linf=[0.41441386613429104])

    analysis_path = joinpath(analysis_directory, "analysis.dat")
    history = accepted_step_history(analysis_path)
    arrays = (; times = history.times, dts = history.dts)
    sources = (analysis_path, [analysis_path, analysis_path], arrays, (history, arrays))
    for (i, source) in enumerate(sources)
        output_path = joinpath(analysis_directory, "time_steps_$i.png")
        fig = plot_time_steps(source; output_path)
        @test fig isa CairoMakie.Figure
        @test isfile(output_path)
        ax = fig.content[1]
        @test ax.xscale[] === identity
        @test ax.yscale[] === identity
        @test length(ax.scene.plots) == (i in (2, 4) ? 2 : 1)
        points = ax.scene.plots[1][1][]
        @test first.(points)≈history.times rtol=1.0e-6
        @test last.(points)≈history.dts rtol=1.0e-6
    end

    output_path = joinpath(analysis_directory, "time_steps_log.pdf")
    fig = plot_time_steps((history, arrays); output_path,
                           labels = ["File history", "Array history"],
                           colors = [:red, :blue], linestyles = [:solid, :dash],
                           legend_position = (:right, :bottom), yscale = log10,
                           xticks = [0.0, 0.5, 1.0], xlims = (0.0, 1.0),
                           yticks = [1.0e-4, 1.0e-2], ylims = (1.0e-5, 1.0))
    @test fig isa CairoMakie.Figure
    @test isfile(output_path)
    @test fig.content[1].yscale[] === log10
end

@testset "convergence triangles on physical spacing" begin
    x = [0.25, 0.5, 1.0]
    series = ((; x, errors = (x .^ 5,), labels = (L"$L^2$",)),)
    mktempdir() do tmpdir
        output_path = joinpath(tmpdir, "spacing_convergence.pdf")
        fig = plot_convergence_1d(series; output_path, xticks = x,
                                  triangle_order = 5, triangle_slope = :positive)
        @test fig isa CairoMakie.Figure
        @test isfile(output_path)
        @test_throws ArgumentError HydroTrixi.plot_bottom_triangle!(fig.content[1], 0.25, 0.5,
                                                                    1.0, 5;
                                                                    triangle_slope = :positive)

        # Non-monotone input order still selects the finest interval. Check a plateau
        # and convergence steeper than the reference slope against both curve endpoints.
        for slope in (:positive, :negative), position in (:below, :above),
            endpoint_errors in ((1.0, 1.0), (1.0e-3, 1.0e-6))
            coarse_error, fine_error = endpoint_errors
            x_values = slope === :positive ? [0.25, 1.0, 0.5] : [4.0, 1.0, 2.0]
            groups = ((; x = x_values,
                       errors = ([fine_error, 1.0, coarse_error],
                                 2 .* [fine_error, 1.0, coarse_error]),
                       labels = (L"$L^2$", L"$L^\infty$")),)
            path = joinpath(tmpdir, "triangle_$(slope)_$(position)_$(fine_error).png")
            fig = plot_convergence_1d(groups; output_path = path, xticks = nothing,
                                      triangle_order = 3, triangle_slope = slope,
                                      triangle_position = position, triangle_gap_factor = 2.0)
            points = fig.content[1].scene.plots[3][1][]
            coarse_x, fine_x = x_values[3], x_values[1]
            @test sort(unique(first.(points))) ≈ sort([coarse_x, fine_x])
            coarse_heights = last.(filter(p -> first(p) == coarse_x, points))
            fine_heights = last.(filter(p -> first(p) == fine_x, points))
            if position === :below
                @test maximum(coarse_heights) <= coarse_error / 2 * (1 + 1.0e-6)
                @test maximum(fine_heights) <= fine_error / 2 * (1 + 1.0e-6)
                @test max(maximum(coarse_heights) / coarse_error,
                          maximum(fine_heights) / fine_error) ≈ 0.5 rtol=1.0e-6
            else
                @test minimum(coarse_heights) >= 4 * coarse_error * (1 - 1.0e-6)
                @test minimum(fine_heights) >= 4 * fine_error * (1 - 1.0e-6)
                @test min(minimum(coarse_heights) / (2 * coarse_error),
                          minimum(fine_heights) / (2 * fine_error)) ≈ 2.0 rtol=1.0e-6
            end
        end
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
