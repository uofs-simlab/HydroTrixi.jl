@inline _scalar_value(u::Number) = u
@inline _scalar_value(u) = u[1]

function plot_solution_1d(sol;
                          output_path = joinpath(pwd(), "solution_1d.pdf"),
                          exact_solution = nothing,
                          numerical_label = LaTeXString("Numerical"),
                          exact_label = LaTeXString("Exact"),
                          xlabel = L"$x$",
                          ylabel = L"$u(x,t)$",
                          title = "",
                          font = DEFAULT_PLOT_FONT,
                          size = DEFAULT_CONVERGENCE_FIGSIZE,
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
                          ylims = nothing,)
    set_serif_tex_theme!(font = font)

    pd = Trixi.PlotData1D(sol; solution_variables = Trixi.cons2cons)
    x = collect(pd.x)
    y = vec(pd.data[:, 1])

    fig = Figure(size = size, fontsize = fontsize)
    ax = Axis(fig[1, 1];
              xlabel = xlabel,
              ylabel = ylabel,
              xlabelfont = xlabelfont,
              ylabelfont = ylabelfont,
              titlefont = titlefont,
              xticklabelfont = xticklabelfont,
              yticklabelfont = yticklabelfont,
              title = isnothing(title) ? "" : title,)

    if !isnothing(xlims)
        CairoMakie.xlims!(ax, xlims)
    end
    if !isnothing(ylims)
        CairoMakie.ylims!(ax, ylims)
    end

    scatterlines!(ax,
                  x,
                  y;
                  label = numerical_label,
                  linewidth = linewidth,
                  markersize = markersize,
                  color = Makie.wong_colors()[1],)

    if !isnothing(exact_solution)
        finite_x = x[isfinite.(x)]
        t_final = sol.t[end]
        x_exact = range(minimum(finite_x), maximum(finite_x); length = 1500)
        y_exact = [_scalar_value(exact_solution(Trixi.SVector(xi), t_final))
                   for xi in x_exact]
        lines!(ax,
               x_exact,
               y_exact;
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
    save(output_path, fig; px_per_unit = 1)

    return fig
end
