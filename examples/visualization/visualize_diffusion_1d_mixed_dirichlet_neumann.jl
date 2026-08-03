# Plot and animate the mixed Dirichlet-Neumann diffusion solution

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi

# Visualization options
final_time = 1.0

Trixi.trixi_include(@__MODULE__,
                    joinpath(dirname(@__DIR__), "elixirs",
                             "elixir_diffusion_1d_mixed_dirichlet_neumann.jl");
                    tspan = (0.0, final_time),
                    saveat = range(0.0, final_time; length = 181))

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
plot_path = joinpath(plots_dir, "diffusion_1d_mixed_solution.pdf")
animation_path = joinpath(plots_dir, "diffusion_1d_mixed_solution.mp4")

plot_solution_1d(sol; exact_solution = exact_solution, output_path = plot_path)
animate_solution_1d(sol; exact_solution = exact_solution, output_path = animation_path,
                    ylims = (0.5, 1.5), framerate = 30)

println("Saved solution plot to: $(plot_path)")
println("Saved solution animation to: $(animation_path)")
