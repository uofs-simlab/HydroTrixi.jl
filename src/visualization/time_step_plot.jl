function HydroTrixi.plot_time_steps(source::Union{AbstractString, NamedTuple, Tuple,
                                                AbstractVector};
                                    output_path = joinpath(pwd(), "time_steps.pdf"),
                                    step_column = "timestep", time_column = "time",
                                    dt_column = "dt", labels = nothing, colors = nothing,
                                    linestyles = nothing, legend_position = (:left, :top),
                                    xlabel = L"$t$ (s)", ylabel = L"$\Delta t$ (s)",
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

    HydroTrixi.set_serif_tex_theme!()
    fig = Figure(size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE, fontsize = 15)
    ax = solution_axis(fig; xlabel, ylabel, xscale, yscale, xticks, yticks, xlims, ylims)
    palette = Makie.wong_colors()
    for (i, history) in enumerate(histories)
        label = isnothing(labels) ? nothing : labels[i]
        color = isnothing(colors) ? palette[mod1(i, length(palette))] : colors[i]
        linestyle = isnothing(linestyles) ? :solid : linestyles[i]
        plot_series!(ax, history.times, history.dts; label, color, linestyle)
    end
    add_legend!(ax; position = legend_position, font = HydroTrixi.DEFAULT_PLOT_FONT,
                labelsize = 14, show_legend = !isnothing(labels))

    mkpath(dirname(abspath(output_path)))
    save(output_path, fig; px_per_unit = 1)
    return fig
end
