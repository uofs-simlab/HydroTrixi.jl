# Example: save the BR1-vs-LDG convergence comparison for an existing mixed
# Dirichlet-Neumann diffusion refinement study.

using CairoMakie
using HydroTrixi
using LaTeXStrings

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))

ldg = results["LDG"]
br1 = results["BR1"]

HydroTrixi.set_serif_tex_theme!()
plot_font = HydroTrixi.DEFAULT_PLOT_FONT

fig = Figure(size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
xticks = HydroTrixi.doubling_dof_ticks(ldg.ndofs; base = 10)
ax = Axis(fig[1, 1];
          xlabel = LaTeXString("Degrees of freedom"),
          ylabel = LaTeXString("Error"),
          xlabelfont = plot_font,
          ylabelfont = plot_font,
          xticklabelfont = plot_font,
          yticklabelfont = plot_font,
          xscale = log10,
          yscale = log10,
          xticks = xticks,)

ax.xminorgridvisible = false
ax.xminorticksvisible = false

colors = Makie.wong_colors()
scatterlines!(ax,
              ldg.ndofs,
              ldg.l2_errors;
              label = L"\mathrm{LDG}\ L^2",
              color = colors[1],
              linestyle = :solid,
              linewidth = 1.8,
              markersize = 7.0)
scatterlines!(ax,
              br1.ndofs,
              br1.l2_errors;
              label = L"\mathrm{BR1}\ L^2",
              color = colors[2],
              linestyle = :solid,
              linewidth = 1.8,
              markersize = 7.0)
scatterlines!(ax,
              ldg.ndofs,
              ldg.linf_errors;
              label = L"\mathrm{LDG}\ L^\infty",
              color = colors[1],
              linestyle = :dash,
              linewidth = 1.8,
              markersize = 7.0)
scatterlines!(ax,
              br1.ndofs,
              br1.linf_errors;
              label = L"\mathrm{BR1}\ L^\infty",
              color = colors[2],
              linestyle = :dash,
              linewidth = 1.8,
              markersize = 7.0)

triangle_order = 4
x_ref_left = ldg.ndofs[end - 1]
x_ref_right = ldg.ndofs[end]
y_ref = min(ldg.l2_errors[end - 1],
            ldg.linf_errors[end - 1],
            br1.l2_errors[end - 1],
            br1.linf_errors[end - 1])
HydroTrixi.plot_bottom_triangle!(ax,
                                 x_ref_left,
                                 x_ref_right,
                                 y_ref,
                                 triangle_order;
                                 triangle_shift = 1.5,
                                 trianglefontsize = 15,
                                 font = plot_font)

axislegend(ax; position = (:right, :top), labelsize = 14, font = plot_font)

comparison_plot = joinpath(plots_dir,
                           "diffusion_equation_1d_mixed_dirichlet_neumann_convergence_br1_vs_ldg.pdf")
save(comparison_plot, fig; px_per_unit = 1)
comparison_plot_png = joinpath(plots_dir,
                               "diffusion_equation_1d_mixed_dirichlet_neumann_convergence_br1_vs_ldg.png")
save(comparison_plot_png, fig; px_per_unit = 1)
println("Saved comparison plot to: $(comparison_plot)")
println("Saved comparison plot to: $(comparison_plot_png)")
