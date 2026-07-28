# Example: save the LDG convergence plot for an existing mixed
# Dirichlet-Neumann diffusion refinement study.

using CairoMakie
using HydroTrixi
using LaTeXStrings

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))

plot_convergence_1d(ndofs, l2_errors, linf_errors;
                    output_path = joinpath(plots_dir,
                                           "diffusion_equation_1d_mixed_dirichlet_neumann_convergence.pdf"),
                    triangle_order = polydeg + 1)
