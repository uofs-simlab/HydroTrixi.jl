# Render an animation of the Celia (1990) Richards-equation infiltration
# profile. Re-runs the elixir with a denser `saveat` so the solution carries
# enough frames for a smooth animation, then animates the pressure-head block
# of the constitutive ODE state (component = 2).

using DiffuSEM
using LaTeXStrings
using Trixi: trixi_include

trixi_include(joinpath(@__DIR__, "elixir_richards_celia_1990.jl");
              saveat = 0.0:2.0:360.0)

plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
if !@isdefined(animation_format)
    animation_format = "mp4"
end

output_path = joinpath(plots_dir,
                       "richards_celia_1990_pressure_head.$(animation_format)")

animate_solution_1d(sol;
                    component = 2,
                    xlabel = L"$z$ (m)",
                    ylabel = L"$\psi$ (m)",
                    ylims = (-0.65, -0.15),
                    output_path = output_path,
                    framerate = 30)

println("Saved Celia ψ(z,t) animation to: $(output_path)")
@show sol.t[end]
