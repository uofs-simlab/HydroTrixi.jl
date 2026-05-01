# Example: save a single final-time comparison plot for an existing mixed
# Dirichlet-Neumann diffusion solution.

using CairoMakie
using HydroTrixi
using LaTeXStrings
plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
if !@isdefined(plot_format)
    plot_format = "pdf"
end

output_path = joinpath(plots_dir,
                       "diffusion_equation_1d_mixed_dirichlet_neumann_solution.$(plot_format)")

plot_solution_1d(sol;
                 exact_solution = exact_solution,
                 output_path = output_path)

println("Saved solution plot to: $(output_path)")
@show sol.t[end]
