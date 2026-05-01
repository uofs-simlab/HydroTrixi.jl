function HydroTrixi.set_serif_tex_theme!(; font = DEFAULT_PLOT_FONT)
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

function apply_axis_limits!(ax; xlims = nothing, ylims = nothing)
    if !isnothing(xlims)
        CairoMakie.xlims!(ax, xlims)
    end
    if !isnothing(ylims)
        CairoMakie.ylims!(ax, ylims)
    end

    return nothing
end
