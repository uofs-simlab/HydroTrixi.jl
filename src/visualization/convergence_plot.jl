function HydroTrixi.doubling_dof_ticks(values::AbstractVector{<:Real}; base::Integer = 10)
    if minimum(values) <= 0
        throw(ArgumentError("Logarithmic DOF axis requires positive values."))
    end
    if base <= 0
        throw(ArgumentError("`base` must be > 0."))
    end

    max_value = maximum(values)
    ticks = Float64[]
    tick = float(base)
    while tick <= max_value * 1.001
        push!(ticks, tick)
        tick *= 2
    end
    isempty(ticks) && push!(ticks, float(base))

    labels = string.(round.(Int, ticks))

    return ticks, labels
end

function HydroTrixi.plot_bottom_triangle!(ax, coarse_x, fine_x, reference_error, order;
                                          gap_factor = 1.5, trianglefontsize = 12,
                                          font = HydroTrixi.DEFAULT_PLOT_FONT,)
    if fine_x <= coarse_x
        throw(ArgumentError("`fine_x` must be greater than `coarse_x`."))
    end
    refinement_ratio = fine_x / coarse_x
    upper_error = reference_error / gap_factor
    lower_error = upper_error / refinement_ratio^order

    # Draw a right triangle with right angle at the bottom-left corner.
    lines!(ax, [coarse_x, fine_x, coarse_x, coarse_x],
           [lower_error, lower_error, upper_error, lower_error]; color = :black,)

    label_x = 10^((2 * log10(coarse_x) + log10(fine_x)) / 3)
    label_error = 10^((2 * log10(lower_error) + log10(upper_error)) / 3)
    text!(ax, label_x, label_error; text = string(order, ":1"),
          align = (:center, :center),
          color = :black, fontsize = trianglefontsize, font = font,)

    return nothing
end

function convergence_triangle_from_data(series_groups, order)
    # Automatic placement assumes every curve uses the same refinement levels.
    reference_x = first(series_groups).x
    if !all(group -> group.x == reference_x, series_groups)
        throw(ArgumentError("Automatic reference triangles require shared x values."))
    end

    # Span the last refinement interval and place the triangle below the lowest error at
    # the coarser of those two levels.
    reference_error = minimum(errors[end - 1] for group in series_groups
                              for errors in group.errors)
    return (; coarse_x = reference_x[end - 1], fine_x = reference_x[end],
            reference_error, order)
end

function HydroTrixi.plot_convergence_1d(series_groups::Union{Tuple, AbstractVector};
                                        output_path = joinpath(pwd(), "convergence_1d.pdf"),
                                        xlabel = LaTeXString("Degrees of freedom"),
                                        ylabel = LaTeXString("Error"),
                                        font = HydroTrixi.DEFAULT_PLOT_FONT,
                                        size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE,
                                        fontsize = 15, legendfontsize = 14,
                                        trianglefontsize = nothing, linewidth = 1.8,
                                        marker = :circle, markersize = 7.0,
                                        show_nodes = true, show_legend = true,
                                        xlabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        ylabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        xticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        yticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        legendfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                        legend_position = (:right, :top), xlims = nothing,
                                        ylims = nothing, xscale = log10, yscale = log10,
                                        xticks = :doubling, triangle_order = nothing,
                                        triangle_gap_factor = 1.5,)
    HydroTrixi.set_serif_tex_theme!(font = font)
    trianglefontsize = isnothing(trianglefontsize) ? fontsize : trianglefontsize

    if xticks === :doubling
        # Cover the complete DOF range when groups contain different refinement levels.
        x_values = reduce(vcat, (collect(group.x) for group in series_groups))
        xticks = HydroTrixi.doubling_dof_ticks(x_values;
                                               base = Int(minimum(x_values)))
    end

    fig = Figure(size = size, fontsize = fontsize)
    ax = solution_axis(fig; xlabel = xlabel, ylabel = ylabel, xlabelfont = xlabelfont,
                       ylabelfont = ylabelfont, xticklabelfont = xticklabelfont,
                       yticklabelfont = yticklabelfont, xscale = xscale,
                       yscale = yscale, xticks = xticks, xlims = xlims, ylims = ylims)

    ax.xminorgridvisible = false
    ax.xminorticksvisible = false

    colors = Makie.wong_colors()
    linestyles = (:solid, :dash, :dot, :dashdot)
    # Share a colour within each formulation and distinguish its error metrics by style
    for (group_index, group) in pairs(series_groups)
        color = get(group, :color, group_index)
        color = color isa Integer ? colors[mod1(color, length(colors))] : color
        for (series_index, errors) in pairs(group.errors)
            linestyle = linestyles[mod1(series_index, length(linestyles))]
            plot_series!(ax, group.x, errors; label = group.labels[series_index],
                         linestyle, linewidth,
                         marker = get(group, :marker, marker),
                         markersize = get(group, :markersize, markersize), color,
                         show_nodes)
        end
    end

    if !isnothing(triangle_order)
        triangle = convergence_triangle_from_data(series_groups, triangle_order)
        HydroTrixi.plot_bottom_triangle!(ax, triangle.coarse_x, triangle.fine_x,
                                         triangle.reference_error, triangle.order;
                                         gap_factor = triangle_gap_factor,
                                         trianglefontsize = trianglefontsize, font = font)
    end

    add_legend!(ax; position = legend_position, font = legendfont,
                labelsize = legendfontsize, show_legend = show_legend)

    outdir = dirname(output_path)
    outdir == "" || mkpath(outdir)
    save(output_path, fig; px_per_unit = 1)

    return fig
end
