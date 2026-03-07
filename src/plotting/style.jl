const DEFAULT_PLOT_FONT = "CMU Serif"

function set_serif_tex_theme!(; font=DEFAULT_PLOT_FONT)
    CairoMakie.set_theme!(CairoMakie.Theme(font=font))
    return nothing
end
