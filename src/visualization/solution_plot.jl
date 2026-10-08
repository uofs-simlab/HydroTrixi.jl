@inline scalar_value(u::Number) = u
@inline scalar_value(u) = u[1]

function exact_solution_x(x)
    finite_x = x[isfinite.(x)]
    return range(minimum(finite_x), maximum(finite_x); length = 1500)
end

function exact_solution_values(exact_solution, x_exact, t)
    return [scalar_value(exact_solution(Trixi.SVector(xi), t)) for xi in x_exact]
end

@inline series_label(show_label, label) = show_label ? label : nothing

function solution_axis(fig; xlabel, ylabel, xlabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                       ylabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                       titlefont = HydroTrixi.DEFAULT_PLOT_FONT,
                       xticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                       yticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT, xscale = identity,
                       yscale = identity, xticks = nothing, yticks = nothing,
                       xlims = nothing,
                       ylims = nothing,)
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, xlabelfont = xlabelfont,
              ylabelfont = ylabelfont, titlefont = titlefont,
              xticklabelfont = xticklabelfont, yticklabelfont = yticklabelfont,
              xscale = xscale, yscale = yscale,)
    if !isnothing(xticks)
        ax.xticks = xticks
    end
    if !isnothing(yticks)
        ax.yticks = yticks
    end
    apply_axis_limits!(ax; xlims = xlims, ylims = ylims)

    return ax
end

function plot_series!(ax, x, y; label = nothing, linestyle = :solid, linewidth = 2.0,
                      marker = :circle, markersize = 7.0, color, show_nodes = false)
    if show_nodes
        scatterlines!(ax, x, y; label = label, linestyle = linestyle, linewidth = linewidth,
                      marker = marker, markersize = markersize, color = color)
    else
        lines!(ax, x, y; label = label, linestyle = linestyle, linewidth = linewidth,
               color = color)
    end

    return nothing
end

function add_legend!(ax; position, font, labelsize, show_legend)
    if show_legend
        axislegend(ax; position = position, font = font, labelsize = labelsize)
    end

    return nothing
end

# Store each curve as points so adaptive meshes update atomically
function solution_points_1d(x, y)
    return Makie.Point2f.(x, y)
end

# Create the shared Makie figure for one-dimensional solution plots
function initialize_solution_plot_1d(x, y, mesh_vertices_x, t;
                                     exact_solution = nothing,
                                     numerical_label = LaTeXString("Numerical"),
                                     exact_label = LaTeXString("Exact"), xlabel = L"$x$",
                                     ylabel = L"$u(x,t)$",
                                     font = HydroTrixi.DEFAULT_PLOT_FONT,
                                     size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE,
                                     fontsize = 15, legendfontsize = 14, linewidth = 2.0,
                                     markersize = 7.0, show_nodes = false,
                                     show_element_boundaries = false,
                                     xlabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                     ylabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                     titlefont = HydroTrixi.DEFAULT_PLOT_FONT,
                                     xticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                     yticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                     legendfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                     legend_position = :rb, xlims = nothing,
                                     ylims = nothing)
    HydroTrixi.set_serif_tex_theme!(font = font)

    show_exact = !isnothing(exact_solution)
    points_obs = Observable(solution_points_1d(x, y))
    mesh_vertices_x_obs = Observable(mesh_vertices_x)
    if isnothing(xlims)
        if isempty(mesh_vertices_x)
            xlims = extrema(filter(isfinite, x))
        else
            xlims = extrema(mesh_vertices_x)
        end
    end

    fig = Figure(size = size, fontsize = fontsize)
    ax = solution_axis(fig; xlabel = xlabel, ylabel = ylabel, xlabelfont = xlabelfont,
                       ylabelfont = ylabelfont, titlefont = titlefont,
                       xticklabelfont = xticklabelfont, yticklabelfont = yticklabelfont,
                       xlims = xlims, ylims = ylims)
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xminorgridvisible = false
    ax.yminorgridvisible = false

    if show_element_boundaries
        vlines!(ax, mesh_vertices_x_obs; color = (:gray, 0.45), linewidth = 0.9)
    end

    if show_nodes
        scatterlines!(ax, points_obs; label = series_label(show_exact, numerical_label),
                      linewidth = linewidth, markersize = markersize,
                      color = Makie.wong_colors()[1])
    else
        lines!(ax, points_obs; label = series_label(show_exact, numerical_label),
               linewidth = linewidth, color = Makie.wong_colors()[1])
    end

    x_exact = nothing
    y_exact_obs = nothing
    if show_exact
        x_exact = exact_solution_x(x)
        y_exact_obs = Observable(exact_solution_values(exact_solution, x_exact, t))
        lines!(ax, x_exact, y_exact_obs; label = exact_label, linewidth = linewidth,
               linestyle = :dash, color = Makie.wong_colors()[2],)
    end
    add_legend!(ax; position = legend_position, font = legendfont,
                labelsize = legendfontsize, show_legend = show_exact)

    return (; fig, points_obs, mesh_vertices_x_obs, x_exact, y_exact_obs, exact_solution)
end

# Update the shared figure for the next saved state
function update_solution_plot_1d!(solution_plot, x, y, mesh_vertices_x, t)
    solution_plot.points_obs[] = solution_points_1d(x, y)
    solution_plot.mesh_vertices_x_obs[] = mesh_vertices_x

    if !isnothing(solution_plot.y_exact_obs)
        solution_plot.y_exact_obs[] = exact_solution_values(solution_plot.exact_solution,
                                                            solution_plot.x_exact, t)
    end

    return nothing
end

function HydroTrixi.plot_solution_1d(data::NamedTuple;
                                     output_path = joinpath(pwd(), "solution_1d.pdf"),
                                     kwargs...)
    solution_plot = initialize_solution_plot_1d(data.x, data.values, data.mesh_vertices_x,
                                                data.time; kwargs...)
    mkpath(dirname(abspath(output_path)))
    save(output_path, solution_plot.fig; px_per_unit = 1)
    return solution_plot.fig
end

function HydroTrixi.plot_solution_1d(sol::SciMLBase.AbstractODESolution;
                                     output_path = joinpath(pwd(), "solution_1d.pdf"),
                                     component = 1, index = lastindex(sol.u),
                                     mesh_history = nothing, kwargs...)
    if isempty(sol.u)
        throw(ArgumentError("The solution has no saved states to plot."))
    end
    if !checkbounds(Bool, sol.u, index)
        throw(BoundsError(sol.u, index))
    end
    if isnothing(mesh_history) && index != lastindex(sol.u)
        throw(ArgumentError("Plotting an earlier saved state requires mesh history."))
    end
    data = HydroTrixi.solution_data_1d(sol; index, component, mesh_history)
    return HydroTrixi.plot_solution_1d(data; output_path, kwargs...)
end

function HydroTrixi.plot_solution_1d(result::HydroTrixi.ImplicitSolveResult;
                                     kwargs...)
    return HydroTrixi.plot_solution_1d(result.sol;
                                       mesh_history = result.mesh_history, kwargs...)
end
