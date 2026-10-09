# Animate the fixed-mesh New Mexico infiltration problem (Celia et al., 1990)

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Trixi

# Simulation and visualization options
final_time = 86_400.0

Trixi.trixi_include(@__MODULE__,
                    joinpath(dirname(@__DIR__), "elixirs",
                             "elixir_richards_celia_new_mexico.jl");
                    amr = false,
                    tspan = (0.0, final_time),
                    saveat = range(0.0, final_time; length = 181))

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
animation_path = joinpath(plots_dir, "richards_celia_new_mexico_pressure_head.mp4")
mass_bias_path = joinpath(plots_dir, "richards_celia_new_mexico_mass_bias.pdf")

animate_solution_1d(sol; component = 2, xlabel = "Distance below surface (m)",
                    ylabel = "Pressure head (m)",
                    ylims = (-10.5, -0.5), output_path = animation_path, framerate = 30)

time_ticks = collect(0.0:20_000.0:80_000.0)
plot_mass_bias_magnitude(sol; output_path = mass_bias_path, yscale = identity,
                         xticks = (time_ticks, string.(Int.(time_ticks))),
                         xlims = (0.0, final_time))

println("Saved New Mexico pressure-head animation to: $(animation_path)")
println("Saved New Mexico mass-bias plot to: $(mass_bias_path)")
