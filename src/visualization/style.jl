const DEFAULT_PLOT_FONT = "CMU Serif"
const DEFAULT_SOLUTION_FIGSIZE = (700, 350)
const DEFAULT_CONVERGENCE_FIGSIZE = (525, 350)

function set_serif_tex_theme!(; font = DEFAULT_PLOT_FONT)
    CairoMakie.set_theme!(CairoMakie.Theme(font = font,
                                           Axis = (xlabelfont = font,
                                                   ylabelfont = font,
                                                   titlefont = font,
                                                   subtitlefont = font,
                                                   xticklabelfont = font,
                                                   yticklabelfont = font),
                                           Legend = (labelfont = font, titlefont = font)))
    return nothing
end

function _apply_axis_limits!(ax; xlims = nothing, ylims = nothing)
    if !isnothing(xlims)
        CairoMakie.xlims!(ax, xlims)
    end
    if !isnothing(ylims)
        CairoMakie.ylims!(ax, ylims)
    end

    return nothing
end
