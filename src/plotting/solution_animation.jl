function animate_solution_1d(sol;
                             output_path = joinpath(pwd(), "solution_1d.mp4"),
                             exact_solution = nothing,
                             numerical_label = LaTeXString("Numerical"),
                             exact_label = LaTeXString("Exact"),
                             xlabel = L"$x$",
                             ylabel = L"$u(x,t)$",
                             show_time = true,
                             font = DEFAULT_PLOT_FONT,
                             size = DEFAULT_SOLUTION_FIGSIZE,
                             fontsize = 15,
                             legendfontsize = 14,
                             linewidth = 2.0,
                             markersize = 7.0,
                             xlabelfont = DEFAULT_PLOT_FONT,
                             ylabelfont = DEFAULT_PLOT_FONT,
                             titlefont = DEFAULT_PLOT_FONT,
                             xticklabelfont = DEFAULT_PLOT_FONT,
                             yticklabelfont = DEFAULT_PLOT_FONT,
                             legendfont = DEFAULT_PLOT_FONT,
                             legend_position = :rb,
                             xlims = nothing,
                             ylims = nothing,
                             framerate = 24,
                             frame_indices = nothing,)
    set_serif_tex_theme!(font = font)

    indices = if isnothing(frame_indices)
        collect(eachindex(sol.t))
    else
        collect(frame_indices)
    end
    isempty(indices) && throw(ArgumentError("`frame_indices` must not be empty."))
    minimum(indices) >= firstindex(sol.t) ||
        throw(ArgumentError("Frame index out of bounds."))
    maximum(indices) <= lastindex(sol.t) ||
        throw(ArgumentError("Frame index out of bounds."))

    semi = sol.prob.p
    first_idx = first(indices)
    pd = Trixi.PlotData1D(sol.u[first_idx], semi; solution_variables = Trixi.cons2cons)
    x = collect(pd.x)
    y = vec(pd.data[:, 1])

    y_obs = Observable(y)
    t_obs = Observable(float(sol.t[first_idx]))

    fig = Figure(size = size, fontsize = fontsize)
    title_obs = if show_time
        @lift "t = $(round($t_obs, sigdigits=5))"
    else
        ""
    end

    ax = Axis(fig[1, 1];
              xlabel = xlabel,
              ylabel = ylabel,
              xlabelfont = xlabelfont,
              ylabelfont = ylabelfont,
              titlefont = titlefont,
              xticklabelfont = xticklabelfont,
              yticklabelfont = yticklabelfont,
              title = title_obs,)

    if !isnothing(xlims)
        CairoMakie.xlims!(ax, xlims)
    end
    if !isnothing(ylims)
        CairoMakie.ylims!(ax, ylims)
    end

    scatterlines!(ax,
                  x,
                  y_obs;
                  label = numerical_label,
                  linewidth = linewidth,
                  markersize = markersize,
                  color = Makie.wong_colors()[1],)

    x_exact = nothing
    y_exact_obs = nothing
    if !isnothing(exact_solution)
        finite_x = x[isfinite.(x)]
        x_exact = range(minimum(finite_x), maximum(finite_x); length = 1500)
        y_exact_obs = Observable([_scalar_value(exact_solution(Trixi.SVector(xi),
                                                               sol.t[first_idx]))
                                  for
                                  xi in x_exact])
        lines!(ax,
               x_exact,
               y_exact_obs;
               label = exact_label,
               linewidth = linewidth,
               linestyle = :dash,
               color = Makie.wong_colors()[2],)
    end

    axislegend(ax;
               position = legend_position,
               font = legendfont,
               labelsize = legendfontsize,)

    outdir = dirname(output_path)
    outdir == "" || mkpath(outdir)

    record(fig, output_path, indices; framerate = framerate) do i
        frame_pd = Trixi.PlotData1D(sol.u[i], semi; solution_variables = Trixi.cons2cons)
        y_obs[] = vec(frame_pd.data[:, 1])
        t_obs[] = float(sol.t[i])

        if !isnothing(y_exact_obs)
            y_exact_obs[] = [_scalar_value(exact_solution(Trixi.SVector(xi), sol.t[i]))
                             for xi in x_exact]
        end
    end

    return output_path
end
