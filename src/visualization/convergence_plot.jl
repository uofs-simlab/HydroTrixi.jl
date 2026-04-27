function compute_eoc(errors::AbstractVector{<:Real}; refinement_factor = 2.0)
    values = similar(errors, Float64)
    fill!(values, NaN)

    denom = log(refinement_factor)

    first_i = firstindex(errors)
    last_i = lastindex(errors)
    for i in (first_i + 1):last_i
        values[i] = log(errors[i - 1] / errors[i]) / denom
    end

    return values
end

function _doubling_dof_ticks(values::AbstractVector{<:Real}; base::Integer = 10)
    minimum(values) > 0 ||
        throw(ArgumentError("Logarithmic DOF axis requires positive values."))
    base > 0 || throw(ArgumentError("`base` must be > 0."))

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

function _plot_bottom_triangle!(ax,
                                x_left,
                                x_right,
                                y_ref,
                                order;
                                triangle_shift = 1.5,
                                trianglefontsize = 12,
                                font = DEFAULT_PLOT_FONT,)
    x_right > x_left || throw(ArgumentError("`x_right` must be greater than `x_left`."))
    ratio = x_right / x_left
    y_top = y_ref / triangle_shift
    y_bottom = y_top / ratio^order

    # Draw a right triangle with right angle at the bottom-left corner.
    lines!(ax,
           [x_left, x_right, x_left, x_left],
           [y_bottom, y_bottom, y_top, y_bottom];
           color = :black,)

    xm = 10^((2 * log10(x_left) + log10(x_right)) / 3)
    ym = 10^((2 * log10(y_bottom) + log10(y_top)) / 3)
    text!(ax,
          xm,
          ym;
          text = string(order, ":1"),
          align = (:center, :center),
          color = :black,
          fontsize = trianglefontsize,
          font = font,)

    return nothing
end

function plot_convergence_1d(ndofs::AbstractVector{<:Real},
                             l2_errors::AbstractVector{<:Real},
                             linf_errors::AbstractVector{<:Real};
                             output_path = joinpath(pwd(), "convergence_1d.pdf"),
                             labels = [L"$L^2$", L"$L^\infty$"],
                             styles = [:solid, :dash],
                             colors = [1, 2],
                             xlabel = LaTeXString("Degrees of freedom"),
                             ylabel = LaTeXString("Error"),
                             font = DEFAULT_PLOT_FONT,
                             size = DEFAULT_CONVERGENCE_FIGSIZE,
                             fontsize = 15,
                             legendfontsize = 14,
                             trianglefontsize = nothing,
                             linewidth = 1.8,
                             xlabelfont = DEFAULT_PLOT_FONT,
                             ylabelfont = DEFAULT_PLOT_FONT,
                             xticklabelfont = DEFAULT_PLOT_FONT,
                             yticklabelfont = DEFAULT_PLOT_FONT,
                             legendfont = DEFAULT_PLOT_FONT,
                             legend_position = (:right, :top),
                             xlims = nothing,
                             ylims = nothing,
                             xscale = log10,
                             yscale = log10,
                             dof_tick_base = 10,
                             triangle_order = nothing,
                             triangle_shift = 1.5,)
    length(ndofs) == length(l2_errors) == length(linf_errors) ||
        throw(ArgumentError("`ndofs`, `l2_errors`, and `linf_errors` must have equal length."))
    length(ndofs) >= 2 ||
        throw(ArgumentError("At least two refinement levels are required."))

    set_serif_tex_theme!(font = font)
    trianglefontsize = isnothing(trianglefontsize) ? fontsize : trianglefontsize

    xticks = _doubling_dof_ticks(ndofs; base = dof_tick_base)

    fig = Figure(size = size, fontsize = fontsize)
    ax = Axis(fig[1, 1];
              xlabel = xlabel,
              ylabel = ylabel,
              xlabelfont = xlabelfont,
              ylabelfont = ylabelfont,
              xticklabelfont = xticklabelfont,
              yticklabelfont = yticklabelfont,
              xscale = xscale,
              yscale = yscale,
              xticks = xticks,)

    ax.xminorgridvisible = false
    ax.xminorticksvisible = false

    _apply_axis_limits!(ax; xlims = xlims, ylims = ylims)

    for (errors, label, style, color) in ((l2_errors, labels[1], styles[1], colors[1]),
                                          (linf_errors, labels[2], styles[2], colors[2]))
        scatterlines!(ax,
                      ndofs,
                      errors;
                      label = label,
                      linestyle = style,
                      linewidth = linewidth,
                      color = Makie.wong_colors()[color],)
    end

    l2_eoc = compute_eoc(l2_errors)
    order = if isnothing(triangle_order)
        finite_orders = filter(isfinite, l2_eoc)
        isempty(finite_orders) ? 1 : round(Int, finite_orders[end])
    else
        triangle_order
    end

    x_ref_left = ndofs[end - 1]
    x_ref_right = ndofs[end]
    y_ref = min(l2_errors[end - 1], linf_errors[end - 1])
    _plot_bottom_triangle!(ax,
                           x_ref_left,
                           x_ref_right,
                           y_ref,
                           order;
                           triangle_shift = triangle_shift,
                           trianglefontsize = trianglefontsize,
                           font = font,)

    axislegend(ax;
               position = legend_position,
               font = legendfont,
               labelsize = legendfontsize,)

    outdir = dirname(output_path)
    outdir == "" || mkpath(outdir)
    save(output_path, fig; px_per_unit = 1)

    return fig
end
