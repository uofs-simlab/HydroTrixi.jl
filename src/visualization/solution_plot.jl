@inline _scalar_value(u::Number) = u
@inline _scalar_value(u) = u[1]

@inline function _plot_component_data(pd, component)
    ncomponents = size(pd.data, 2)
    1 <= component <= ncomponents ||
        throw(ArgumentError("Component $component is out of bounds for $ncomponents variables."))
    return vec(pd.data[:, component])
end

function _plot_data_1d(u_ode, semi; component = 1)
    if semi isa SemidiscretizationConstitutive
        nstate_variables = Trixi.nvariables(semi.semidiscretization)
        1 <= component <= 2 * nstate_variables ||
            throw(ArgumentError("Component $component is out of bounds for $(2 * nstate_variables) variables."))
        if component <= nstate_variables
            data = evolved_variable_view(u_ode)
            local_component = component
        else
            data = state_variable_view(u_ode)
            local_component = component - nstate_variables
        end
        pd = Trixi.PlotData1D(data, semi.semidiscretization;
                              solution_variables = Trixi.cons2cons)
        return pd, local_component
    end

    return Trixi.PlotData1D(u_ode, semi; solution_variables = Trixi.cons2cons), component
end

function _plot_data_1d(sol; component = 1)
    semi = sol.prob.p
    if semi isa SemidiscretizationConstitutive
        return _plot_data_1d(sol.u[end], semi; component = component)
    end

    return Trixi.PlotData1D(sol; solution_variables = Trixi.cons2cons), component
end

function _plot_curve_1d(u_ode, semi; component = 1)
    pd, local_component = _plot_data_1d(u_ode, semi; component = component)
    return collect(pd.x), _plot_component_data(pd, local_component)
end

function _plot_curve_1d(sol; component = 1)
    pd, local_component = _plot_data_1d(sol; component = component)
    return collect(pd.x), _plot_component_data(pd, local_component)
end

function _exact_solution_x(x)
    finite_x = x[isfinite.(x)]
    return range(minimum(finite_x), maximum(finite_x); length = 1500)
end

function _exact_solution_values(exact_solution, x_exact, t)
    return [_scalar_value(exact_solution(Trixi.SVector(xi), t)) for xi in x_exact]
end

function _solution_axis(fig;
                        xlabel,
                        ylabel,
                        title = "",
                        xlabelfont = DEFAULT_PLOT_FONT,
                        ylabelfont = DEFAULT_PLOT_FONT,
                        titlefont = DEFAULT_PLOT_FONT,
                        xticklabelfont = DEFAULT_PLOT_FONT,
                        yticklabelfont = DEFAULT_PLOT_FONT,
                        xlims = nothing,
                        ylims = nothing,)
    ax = Axis(fig[1, 1];
              xlabel = xlabel,
              ylabel = ylabel,
              xlabelfont = xlabelfont,
              ylabelfont = ylabelfont,
              titlefont = titlefont,
              xticklabelfont = xticklabelfont,
              yticklabelfont = yticklabelfont,
              title = title,)
    _apply_axis_limits!(ax; xlims = xlims, ylims = ylims)

    return ax
end

function plot_solution_1d(sol;
                          output_path = joinpath(pwd(), "solution_1d.pdf"),
                          component = 1,
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

    x, y = _plot_curve_1d(sol; component = component)

    fig = Figure(size = size, fontsize = fontsize)
    ax = _solution_axis(fig;
                        xlabel = xlabel,
                        ylabel = ylabel,
                        title = title,
                        xlabelfont = xlabelfont,
                        ylabelfont = ylabelfont,
                        titlefont = titlefont,
                        xticklabelfont = xticklabelfont,
                        yticklabelfont = yticklabelfont,
                        xlims = xlims,
                        ylims = ylims)

    scatterlines!(ax,
                  x,
                  y;
                  label = numerical_label,
                  linewidth = linewidth,
                  markersize = markersize,
                  color = Makie.wong_colors()[1],)

    if !isnothing(exact_solution)
        t_final = sol.t[end]
        x_exact = _exact_solution_x(x)
        y_exact = _exact_solution_values(exact_solution, x_exact, t_final)
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
