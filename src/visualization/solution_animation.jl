function animate_solution_1d(sol;
                             output_path = joinpath(pwd(), "solution_1d.mp4"),
                             component = 1,
                             exact_solution = nothing,
                             numerical_label = LaTeXString("Numerical"),
                             exact_label = LaTeXString("Exact"),
                             xlabel = L"$x$",
                             ylabel = L"$u(x,t)$",
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
    x, y = _plot_curve_1d(sol.u[first_idx], semi; component = component)

    y_obs = Observable(y)

    fig = Figure(size = size, fontsize = fontsize)

    ax = _solution_axis(fig;
                        xlabel = xlabel,
                        ylabel = ylabel,
                        title = "",
                        xlabelfont = xlabelfont,
                        ylabelfont = ylabelfont,
                        titlefont = titlefont,
                        xticklabelfont = xticklabelfont,
                        yticklabelfont = yticklabelfont,
                        xlims = xlims,
                        ylims = ylims)

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
        x_exact = _exact_solution_x(x)
        y_exact_obs = Observable(_exact_solution_values(exact_solution,
                                                        x_exact,
                                                        sol.t[first_idx]))
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
        _, y_frame = _plot_curve_1d(sol.u[i], semi; component = component)
        y_obs[] = y_frame

        if !isnothing(y_exact_obs)
            y_exact_obs[] = _exact_solution_values(exact_solution, x_exact, sol.t[i])
        end
    end

    return output_path
end
