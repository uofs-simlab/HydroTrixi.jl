# Tuples and vectors collect multiple solutions or analysis files.
@inline function is_mass_bias_source_collection(source)
    if source isa Tuple
        return true
    end
    if !(source isa AbstractVector)
        return false
    end
    return !(hasproperty(source, :prob) && hasproperty(source, :u))
end

function mass_bias_log_axis(histories)
    positive_values = Float64[]
    for (_, biases) in histories, bias in biases
        magnitude = abs(bias)
        if isfinite(magnitude) && magnitude > 0
            push!(positive_values, magnitude)
        end
    end
    if isempty(positive_values)
        throw(ArgumentError("Mass-bias histories contain no positive finite magnitudes"))
    end

    lower_exponent = floor(Int, log10(minimum(positive_values)))
    upper_exponent = ceil(Int, log10(maximum(positive_values)))
    if lower_exponent == upper_exponent
        upper_exponent += 1
    end
    tick_exponents = filter(iseven, lower_exponent:upper_exponent)
    if length(tick_exponents) < 2
        tick_exponents = collect(lower_exponent:upper_exponent)
    end
    ticks = 10.0 .^ tick_exponents
    tick_labels = [LaTeXString("10^{$exponent}") for exponent in tick_exponents]
    return (; limits = (10.0^lower_exponent, 10.0^upper_exponent),
            ticks = (ticks, tick_labels))
end

function mass_bias_histories(source; time_column, mass_balance_column)
    source_items = is_mass_bias_source_collection(source) ? source : (source,)
    return map(source_items) do source_item
        if source_item isa AbstractString
            HydroTrixi.mass_bias_history(source_item; time_column,
                                          mass_balance_column)
        else
            HydroTrixi.mass_bias_history(source_item)
        end
    end
end

function HydroTrixi.mass_bias_magnitude_axis(source; time_column = "time",
                                             mass_balance_column = "mass_balance")
    histories = mass_bias_histories(source; time_column, mass_balance_column)
    return mass_bias_log_axis(histories)
end

function HydroTrixi.plot_mass_bias_magnitude(source;
                                             output_path = joinpath(pwd(), "mass_bias.pdf"),
                                             time_column = "time",
                                             mass_balance_column = "mass_balance",
                                             labels = nothing, colors = nothing,
                                             linestyles = nothing,
                                             legend_position = (:right, :top),
                                             yscale = log10,
                                             xticks = nothing, yticks = nothing,
                                             xlims = nothing, ylims = nothing)
    if yscale !== log10 && yscale !== identity
        throw(ArgumentError("`yscale` must be `log10` or `identity`"))
    end
    source_items = is_mass_bias_source_collection(source) ? source : (source,)
    for (option, values) in ((:labels, labels), (:colors, colors),
                             (:linestyles, linestyles))
        if !isnothing(values) && length(values) != length(source_items)
            throw(ArgumentError("`$option` must have one entry per source"))
        end
    end

    histories = mass_bias_histories(source_items; time_column, mass_balance_column)
    magnitudes = map(histories) do (_, biases)
        values = abs.(biases)
        if yscale === log10
            return map(value -> value > 0 ? value : NaN, values)
        end
        return values
    end
    if yscale === log10 && (isnothing(ylims) || isnothing(yticks))
        bias_axis = mass_bias_log_axis(histories)
        if isnothing(ylims)
            ylims = bias_axis.limits
        end
        if isnothing(yticks)
            yticks = bias_axis.ticks
        end
    elseif yscale === identity && (isnothing(ylims) || isnothing(yticks))
        finite_values = filter(isfinite, vcat(magnitudes...))
        if isempty(finite_values)
            throw(ArgumentError("Mass-bias histories contain no finite magnitudes"))
        end
        maximum_value = isnothing(ylims) ? maximum(finite_values) : ylims[2]
        if maximum_value > 0
            exponent = floor(Int, log10(maximum_value))
            scale = 10.0^exponent
            normalized_maximum = maximum_value / scale
            tick_step = if normalized_maximum <= 2
                0.5
            elseif normalized_maximum <= 5
                1.0
            else
                2.0
            end
            upper_factor = ceil(normalized_maximum / tick_step) * tick_step
            if isnothing(ylims)
                ylims = (0.0, upper_factor * scale)
            end
            if isnothing(yticks)
                factors = collect(0.0:tick_step:upper_factor)
                tick_labels = map(factors) do factor
                    if iszero(factor)
                        return LaTeXString("0")
                    end
                    coefficient = isinteger(factor) ? string(Int(factor)) : string(factor)
                    return LaTeXString("$(coefficient)\\times 10^{$exponent}")
                end
                yticks = (scale .* factors, tick_labels)
            end
        else
            if isnothing(ylims)
                ylims = (0.0, 1.0)
            end
            if isnothing(yticks)
                yticks = ([0.0, 1.0], ["0", "1"])
            end
        end
    end

    HydroTrixi.set_serif_tex_theme!()
    fig = Figure(size = HydroTrixi.DEFAULT_SOLUTION_FIGSIZE, fontsize = 15)
    ax = solution_axis(fig; xlabel = L"$t$ (s)",
                       ylabel = L"$|\epsilon_{\mathrm{b}}(t)|$ (m)",
                       yscale, xticks, yticks, xlims, ylims)

    palette = Makie.wong_colors()
    for (i, ((times, _), values)) in enumerate(zip(histories, magnitudes))
        label = isnothing(labels) ? nothing : labels[i]
        color = isnothing(colors) ? palette[mod1(i, length(palette))] : colors[i]
        linestyle = isnothing(linestyles) ? :solid : linestyles[i]
        plot_series!(ax, times, values; label, color, linestyle)
    end

    add_legend!(ax; position = legend_position, font = HydroTrixi.DEFAULT_PLOT_FONT,
                labelsize = 14, show_legend = !isnothing(labels))

    output_directory = dirname(output_path)
    if output_directory != ""
        mkpath(output_directory)
    end
    save(output_path, fig; px_per_unit = 1)
    return fig
end
