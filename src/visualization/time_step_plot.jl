function HydroTrixi.plot_time_steps(source::Union{AbstractString, NamedTuple, Tuple,
                                                AbstractVector};
                                    output_path = joinpath(pwd(), "time_steps.pdf"),
                                    step_column = "timestep", time_column = "time",
                                    dt_column = "dt", labels = nothing, colors = nothing,
                                    linestyles = nothing,
                                    font = HydroTrixi.DEFAULT_PLOT_FONT,
                                    size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE,
                                    fontsize = 15, legendfontsize = 14, linewidth = 2.0,
                                    show_legend = !isnothing(labels),
                                    xlabelfont = font, ylabelfont = font,
                                    titlefont = font, xticklabelfont = font,
                                    yticklabelfont = font, legendfont = font,
                                    legend_position = (:left, :top),
                                    xlabel = L"Time $t$ (s)",
                                    ylabel = L"Time step size $\Delta t$ (s)",
                                    xscale = identity, yscale = identity,
                                    xticks = nothing, yticks = nothing,
                                    xlims = nothing, ylims = nothing)
    source_items = if source isa AbstractString || source isa NamedTuple
        (source,)
    else
        source
    end
    for (option, values) in ((:labels, labels), (:colors, colors),
                             (:linestyles, linestyles))
        if !isnothing(values) && length(values) != length(source_items)
            throw(ArgumentError("`$option` must have one entry per source"))
        end
    end

    histories = map(source_items) do source_item
        if source_item isa AbstractString
            return HydroTrixi.accepted_step_history(source_item; step_column, time_column,
                                                    dt_column)
        end
        return source_item
    end

    HydroTrixi.set_serif_tex_theme!(; font)
    fig = Figure(; size, fontsize)
    ax = solution_axis(fig; xlabel, ylabel, xlabelfont, ylabelfont, titlefont,
                       xticklabelfont, yticklabelfont, xscale, yscale, xticks, yticks,
                       xlims, ylims)
    palette = Makie.wong_colors()
    for (i, history) in enumerate(histories)
        label = isnothing(labels) ? nothing : labels[i]
        color = isnothing(colors) ? palette[mod1(i, length(palette))] : colors[i]
        linestyle = isnothing(linestyles) ? :solid : linestyles[i]
        plot_series!(ax, history.times, history.dts; label, color, linestyle, linewidth)
    end
    add_legend!(ax; position = legend_position, font = legendfont,
                labelsize = legendfontsize, show_legend)

    mkpath(dirname(abspath(output_path)))
    save(output_path, fig; px_per_unit = 1)
    return fig
end
