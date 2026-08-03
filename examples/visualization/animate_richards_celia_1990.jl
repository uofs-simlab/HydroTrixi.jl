# Animate the fixed-mesh Celia benchmark and plot its mass bias

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi

# Simulation and visualization options
final_time = 360.0

Trixi.trixi_include(@__MODULE__,
                    joinpath(dirname(@__DIR__), "elixirs",
                             "elixir_richards_celia_1990.jl");
                    tspan = (0.0, final_time),
                    saveat = range(0.0, final_time; length = 181))

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
animation_path = joinpath(plots_dir, "richards_celia_1990_pressure_head.mp4")
mass_bias_path = joinpath(plots_dir, "richards_celia_1990_mass_bias.pdf")

animate_solution_1d(sol; component = 2, xlabel = L"$z$ (m)", ylabel = L"$\psi$ (m)",
                    ylims = (-0.65, -0.15), output_path = animation_path, framerate = 30)

plot_mass_bias(sol; output_path = mass_bias_path, xlabel = L"$t$ (s)",
               ylabel = L"$\epsilon_b$ (m)",
               xticks = range(0.0, final_time; length = 7))

println("Saved Celia pressure-head animation to: $(animation_path)")
println("Saved Celia mass-bias plot to: $(mass_bias_path)")
