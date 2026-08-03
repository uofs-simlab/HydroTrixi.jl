# Compare LDG and BR1 for the mixed Dirichlet-Neumann diffusion problem

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Printf: @printf
using Trixi

elixir = joinpath(dirname(@__DIR__), "elixirs",
                  "elixir_diffusion_1d_mixed_dirichlet_neumann.jl")

# Convergence study options
diffusivity = 0.5
forcing_amplitude = 0.4
forcing_frequency = 4.0
polydeg = 3
tspan = (0.0, 0.25)
base_initial_refinement_level = 2
iterations = 4

schemes = ((; name = "LDG", solver_parabolic = ParabolicFormulationLocalDG(), color = 1),
           (; name = "BR1", solver_parabolic = ParabolicFormulationBassiRebay1(),
            color = 2))

function run_scheme(scheme)
    _, errors = convergence_test(@__MODULE__, elixir, iterations;
                                 diffusivity = diffusivity,
                                 forcing_amplitude = forcing_amplitude,
                                 angular_frequency = 2 * pi * forcing_frequency,
                                 dirichlet_mean = 1.0, polydeg = polydeg,
                                 tspan = tspan,
                                 initial_refinement_level = base_initial_refinement_level,
                                 solver_parabolic = scheme.solver_parabolic)

    levels = collect(base_initial_refinement_level:
                     (base_initial_refinement_level + iterations - 1))
    return (; levels, ndofs = (polydeg + 1) .* 2 .^ levels,
            l2_errors = vec(errors[:l2]), linf_errors = vec(errors[:linf]))
end

results = Dict(scheme.name => run_scheme(scheme) for scheme in schemes)

println("1D diffusion convergence study with mixed Dirichlet-Neumann boundaries")
println("diffusivity = $(diffusivity), polydeg = $(polydeg), tspan = $(tspan)")
println("forcing_amplitude = $(forcing_amplitude), forcing_frequency = $(forcing_frequency)")
println()

for scheme in schemes
    data = results[scheme.name]
    l2_eoc = compute_eoc(data.l2_errors)
    linf_eoc = compute_eoc(data.linf_errors)
    println("Scheme: $(scheme.name)")
    @printf("%-7s %-10s %-14s %-8s %-14s %-8s\n", "level", "ndofs", "L2 error",
            "EOC", "Linf error", "EOC")
    for i in eachindex(data.levels)
        @printf("%-7d %-10d %-14.6e %-8.3f %-14.6e %-8.3f\n", data.levels[i],
                data.ndofs[i], data.l2_errors[i], l2_eoc[i], data.linf_errors[i],
                linf_eoc[i])
    end
    println()
end

ldg = results["LDG"]
br1 = results["BR1"]
println("BR1/LDG error ratio by level")
@printf("%-7s %-12s %-12s\n", "level", "L2 ratio", "Linf ratio")
for i in eachindex(ldg.levels)
    @printf("%-7d %-12.4f %-12.4f\n", ldg.levels[i],
            br1.l2_errors[i]/ldg.l2_errors[i],
            br1.linf_errors[i]/ldg.linf_errors[i])
end

# Plot both parabolic formulations on one axis
HydroTrixi.set_serif_tex_theme!()
plot_font = HydroTrixi.DEFAULT_PLOT_FONT
figure = Figure(size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
axis = Axis(figure[1, 1]; xlabel = LaTeXString("Degrees of freedom"),
            ylabel = LaTeXString("Error"), xlabelfont = plot_font,
            ylabelfont = plot_font, xticklabelfont = plot_font,
            yticklabelfont = plot_font, xscale = log10, yscale = log10,
            xticks = HydroTrixi.doubling_dof_ticks(ldg.ndofs;
                                                   base = minimum(ldg.ndofs)))
axis.xminorgridvisible = false
axis.xminorticksvisible = false

colors = Makie.wong_colors()
for scheme in schemes
    data = results[scheme.name]
    scatterlines!(axis, data.ndofs, data.l2_errors; label = "$(scheme.name) L²",
                  color = colors[scheme.color], linestyle = :solid, linewidth = 1.8)
    scatterlines!(axis, data.ndofs, data.linf_errors; label = "$(scheme.name) L∞",
                  color = colors[scheme.color], linestyle = :dash, linewidth = 1.8)
end

y_ref = min(ldg.l2_errors[end - 1], ldg.linf_errors[end - 1],
            br1.l2_errors[end - 1], br1.linf_errors[end - 1])
HydroTrixi.plot_bottom_triangle!(axis, ldg.ndofs[end - 1], ldg.ndofs[end], y_ref,
                                 polydeg + 1; trianglefontsize = 15, font = plot_font)
axislegend(axis; position = (:right, :top), labelsize = 14, font = plot_font)

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
output_path = joinpath(plots_dir, "diffusion_1d_mixed_br1_vs_ldg.pdf")
save(output_path, figure; px_per_unit = 1)
println("Saved convergence plot to: $(output_path)")
