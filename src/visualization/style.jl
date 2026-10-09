function HydroTrixi.set_serif_tex_theme!(; font = HydroTrixi.DEFAULT_PLOT_FONT)
    # Match ordinary text to the bundled font family used for LaTeX labels
    CairoMakie.set_theme!(CairoMakie.theme_latexfonts(); font = font,
                          Axis = (xlabelfont = font, ylabelfont = font,
                                  titlefont = font, subtitlefont = font,
                                  xticklabelfont = font, yticklabelfont = font),
                          Legend = (labelfont = font, titlefont = font))
    return nothing
end
