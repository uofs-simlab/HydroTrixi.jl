# Tuples and vectors collect multiple solutions or analysis files
@inline function is_mass_bias_source_collection(source)
    source isa Tuple && return true
    source isa AbstractVector || return false
    return !(hasproperty(source, :prob) && hasproperty(source, :u))
end

function HydroTrixi.plot_mass_bias(source; output_path = joinpath(pwd(), "mass_bias.pdf"),
                                   initial_water_content = nothing, time_column = "time",
                                   mass_bias_column = "mass_bias", label = L"$\epsilon_b$",
                                   labels = nothing, xlabel = L"$t$",
                                   ylabel = L"$\epsilon_b$", absolute = false,
                                   ynorm = :auto, exponent_text = :auto,
                                   font = HydroTrixi.DEFAULT_PLOT_FONT,
                                   size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE,
                                   fontsize = 15, legendfontsize = 14, linewidth = 2.0,
                                   markersize = 7.0, show_nodes = false,
                                   show_legend = is_mass_bias_source_collection(source),
                                   xlabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                   ylabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                   titlefont = HydroTrixi.DEFAULT_PLOT_FONT,
                                   xticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                   yticklabelfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                   legendfont = HydroTrixi.DEFAULT_PLOT_FONT,
                                   legend_position = :rb, yscale = identity,
                                   xticks = nothing, xlims = nothing, ylims = nothing,)
    HydroTrixi.set_serif_tex_theme!(font = font)

    # Load each requested mass-bias history from a solution or analysis file
    source_items = is_mass_bias_source_collection(source) ? source : (source,)
    histories = map(source_items) do source_item
        if source_item isa AbstractString
            if !isnothing(initial_water_content)
                throw(ArgumentError("`initial_water_content` is only supported for " *
                                    "saved solutions."))
            end
            HydroTrixi.mass_bias_history(source_item; time_column = time_column,
                                         mass_bias_column = mass_bias_column)
        else
            HydroTrixi.mass_bias_history(source_item;
                                         initial_water_content = initial_water_content)
        end
    end
    times = first.(histories)
    biases = last.(histories)

    # Apply display transforms before optional exponent scaling
    biases = absolute ? [abs.(bias) for bias in biases] : biases
    if yscale === log10
        biases = [map(value -> value > 0 ? value : NaN, bias) for bias in biases]
    end

    # Choose a shared exponent for all finite nonzero mass-bias values
    if ynorm === :auto
        finite_biases = filter(isfinite, abs.(vcat(biases...)))
        filter!(!iszero, finite_biases)
        exponent = isempty(finite_biases) ? 0 :
                   floor(Int, log10(maximum(finite_biases)))
        ynorm = 10.0^exponent
        if exponent_text === :auto
            exponent_text = exponent == 0 ? nothing :
                            LaTeXString("\\times\\!10^{$exponent}")
        end
    elseif exponent_text === :auto
        exponent_text = nothing
    end
    biases = [bias ./ ynorm for bias in biases]

    series_labels = if !isnothing(labels)
        labels
    elseif length(histories) == 1
        (label,)
    else
        [LaTeXString("Run $i") for i in eachindex(histories)]
    end

    fig = Figure(size = size, fontsize = fontsize)
    ax = solution_axis(fig; xlabel = xlabel, ylabel = ylabel, xlabelfont = xlabelfont,
                       ylabelfont = ylabelfont, titlefont = titlefont,
                       xticklabelfont = xticklabelfont, yticklabelfont = yticklabelfont,
                       yscale = yscale, xticks = xticks, xlims = xlims, ylims = ylims)

    if !isnothing(exponent_text)
        Label(fig[1, 1, Top()], halign = :left, exponent_text)
    end

    # Draw each mass-bias history with a deterministic palette entry
    colors = Makie.wong_colors()
    for i in eachindex(histories)
        plot_series!(ax, times[i], biases[i];
                     label = series_label(show_legend, series_labels[i]),
                     linewidth = linewidth, markersize = markersize,
                     color = colors[mod1(i, length(colors))], show_nodes = show_nodes)
    end

    add_legend!(ax; position = legend_position, font = legendfont,
                labelsize = legendfontsize, show_legend = show_legend)

    outdir = dirname(output_path)
    outdir == "" || mkpath(outdir)
    save(output_path, fig; px_per_unit = 1)

    return fig
end
