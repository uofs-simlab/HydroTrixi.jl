# Choose a shared exponent for all finite nonzero mass-bias values
function mass_bias_exponent(values)
    finite_values = filter(isfinite, abs.(values))
    filter!(!iszero, finite_values)
    isempty(finite_values) && return 0
    return floor(Int, log10(maximum(finite_values)))
end

function mass_bias_exponent_text(exponent)
    exponent == 0 && return nothing
    return LaTeXString("\\times\\!10^{$exponent}")
end

function scaled_mass_biases(values, ynorm, exponent_text)
    if ynorm === :auto
        exponent = mass_bias_exponent(vcat(values...))
        ynorm = 10.0^exponent
        scaled_values = [value ./ ynorm for value in values]
        exponent_text === :auto && return scaled_values, mass_bias_exponent_text(exponent)
    end

    scaled_values = [value ./ ynorm for value in values]
    exponent_text === :auto && return scaled_values, nothing
    return scaled_values, exponent_text
end

# Apply display transforms before optional exponent scaling
function transformed_mass_biases(values, absolute, yscale)
    transformed_values = absolute ? [abs.(value) for value in values] : values
    yscale === log10 || return transformed_values
    return [map(value -> value > 0 ? value : NaN, series)
            for series in transformed_values]
end

# SciML solutions are single plot sources, not source collections
@inline function is_mass_bias_saved_solution(source)
    return hasproperty(source, :prob) && hasproperty(source, :u)
end

# Tuples and vectors collect multiple solutions or analysis files
@inline function is_mass_bias_source_collection(source)
    source isa Tuple && return true
    source isa AbstractVector || return false
    return !is_mass_bias_saved_solution(source)
end

# Normalize one source to the iterable shape used for multiple sources
@inline function mass_bias_source_items(source)
    is_mass_bias_source_collection(source) && return source
    return (source,)
end

# Read one mass-bias history from a saved solution or analysis file
function mass_bias_series(source; initial_water_content, time_column, mass_bias_column)
    if source isa AbstractString
        isnothing(initial_water_content) ||
            throw(ArgumentError("`initial_water_content` is only supported for " *
                                "saved solutions."))
        return HydroTrixi.mass_bias_history(source;
                                            time_column = time_column,
                                            mass_bias_column = mass_bias_column)
    end

    return HydroTrixi.mass_bias_history(source;
                                        initial_water_content = initial_water_content)
end

function mass_bias_labels(nseries, label, labels)
    !isnothing(labels) && return labels
    nseries == 1 && return (label,)
    return [LaTeXString("Run $i") for i in 1:nseries]
end

function HydroTrixi.plot_mass_bias(source;
                                   output_path = joinpath(pwd(), "mass_bias.pdf"),
                                   initial_water_content = nothing,
                                   time_column = "time",
                                   mass_bias_column = "mass_bias",
                                   label = L"$\epsilon_b$",
                                   labels = nothing,
                                   xlabel = L"$t$",
                                   ylabel = L"$\epsilon_b$",
                                   absolute = false,
                                   ynorm = :auto,
                                   exponent_text = :auto,
                                   font = DEFAULT_PLOT_FONT,
                                   size = DEFAULT_SOLUTION_FIGSIZE,
                                   fontsize = 15,
                                   legendfontsize = 14,
                                   linewidth = 2.0,
                                   markersize = 7.0,
                                   show_nodes = false,
                                   show_legend = is_mass_bias_source_collection(source),
                                   xlabelfont = DEFAULT_PLOT_FONT,
                                   ylabelfont = DEFAULT_PLOT_FONT,
                                   titlefont = DEFAULT_PLOT_FONT,
                                   xticklabelfont = DEFAULT_PLOT_FONT,
                                   yticklabelfont = DEFAULT_PLOT_FONT,
                                   legendfont = DEFAULT_PLOT_FONT,
                                   legend_position = :rb,
                                   yscale = identity,
                                   xticks = nothing,
                                   xlims = nothing,
                                   ylims = nothing,)
    set_serif_tex_theme!(font = font)

    # Load and transform all requested mass-bias histories
    histories = [mass_bias_series(source_item;
                                  initial_water_content = initial_water_content,
                                  time_column = time_column,
                                  mass_bias_column = mass_bias_column)
                 for source_item in mass_bias_source_items(source)]
    times = first.(histories)
    biases = last.(histories)
    biases = transformed_mass_biases(biases, absolute, yscale)
    biases, exponent_text = scaled_mass_biases(biases, ynorm, exponent_text)
    series_labels = mass_bias_labels(length(histories), label, labels)

    fig = Figure(size = size, fontsize = fontsize)
    ax = solution_axis(fig;
                       xlabel = xlabel,
                       ylabel = ylabel,
                       xlabelfont = xlabelfont,
                       ylabelfont = ylabelfont,
                       titlefont = titlefont,
                       xticklabelfont = xticklabelfont,
                       yticklabelfont = yticklabelfont,
                       yscale = yscale,
                       xticks = xticks,
                       xlims = xlims,
                       ylims = ylims)

    if !isnothing(exponent_text)
        Label(fig[1, 1, Top()], halign = :left, exponent_text)
    end

    # Draw each mass-bias history with a deterministic palette entry
    colors = Makie.wong_colors()
    for i in eachindex(histories)
        plot_series!(ax,
                     times[i],
                     biases[i];
                     label = series_label(show_legend, series_labels[i]),
                     linewidth = linewidth,
                     markersize = markersize,
                     color = colors[mod1(i, length(colors))],
                     show_nodes = show_nodes)
    end

    add_legend!(ax;
                position = legend_position,
                font = legendfont,
                labelsize = legendfontsize,
                show_legend = show_legend)

    outdir = dirname(output_path)
    outdir == "" || mkpath(outdir)
    save(output_path, fig; px_per_unit = 1)

    return fig
end
