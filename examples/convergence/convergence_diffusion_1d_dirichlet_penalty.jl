# Compare ordinary and penalty Dirichlet boundaries for the LDG diffusion discretization

using CairoMakie
using HydroTrixi
using LaTeXStrings
using Printf: @printf
using Trixi

elixir = joinpath(dirname(@__DIR__), "elixirs",
                  "elixir_diffusion_1d_dirichlet_dirichlet.jl")

# Convergence study options
diffusivity = 0.5
polydeg = 3
tspan = (0.0, 0.25)
base_initial_refinement_level = 2
iterations = 5

schemes = ((; name = "LDG", use_boundary_penalty = false, color = 1),
           (; name = "LDG + boundary penalty", use_boundary_penalty = true, color = 2))

function run_scheme(scheme)
    _, errors = convergence_test(@__MODULE__, elixir, iterations;
                                 diffusivity = diffusivity, polydeg = polydeg,
                                 tspan = tspan,
                                 initial_refinement_level = base_initial_refinement_level,
                                 solver_parabolic = ParabolicFormulationLocalDG(),
                                 use_boundary_penalty = scheme.use_boundary_penalty)

    levels = collect(base_initial_refinement_level:
                     (base_initial_refinement_level + iterations - 1))
    return (; levels, ndofs = (polydeg + 1) .* 2 .^ levels,
            dts = 0.01 .* (1.0 ./ (2.0 .^ levels)) .^ 2,
            l2_errors = vec(errors[:l2]), linf_errors = vec(errors[:linf]))
end

results = Dict(scheme.name => run_scheme(scheme) for scheme in schemes)

println("1D diffusion convergence study with double Dirichlet boundaries")
println("diffusivity = $(diffusivity), polydeg = $(polydeg), tspan = $(tspan)")
println()

for scheme in schemes
    data = results[scheme.name]
    l2_eoc = compute_eoc(data.l2_errors)
    linf_eoc = compute_eoc(data.linf_errors)
    println("Scheme: $(scheme.name)")
    @printf("%-7s %-10s %-11s %-14s %-8s %-14s %-8s\n", "level", "ndofs", "dt",
            "L2 error", "EOC", "Linf error", "EOC")
    for i in eachindex(data.levels)
        @printf("%-7d %-10d %-11.3e %-14.6e %-8.3f %-14.6e %-8.3f\n",
                data.levels[i], data.ndofs[i], data.dts[i], data.l2_errors[i],
                l2_eoc[i], data.linf_errors[i], linf_eoc[i])
    end
    println()
end

# Plot both boundary treatments on one axis
HydroTrixi.set_serif_tex_theme!()
plot_font = HydroTrixi.DEFAULT_PLOT_FONT
reference_data = results[first(schemes).name]
figure = Figure(size = HydroTrixi.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
axis = Axis(figure[1, 1]; xlabel = LaTeXString("Degrees of freedom"),
            ylabel = LaTeXString("Error"), xlabelfont = plot_font,
            ylabelfont = plot_font, xticklabelfont = plot_font,
            yticklabelfont = plot_font, xscale = log10, yscale = log10,
            xticks = HydroTrixi.doubling_dof_ticks(reference_data.ndofs;
                                                   base = minimum(reference_data.ndofs)))
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

y_ref = minimum(min(results[scheme.name].l2_errors[end - 1],
                    results[scheme.name].linf_errors[end - 1]) for scheme in schemes)
HydroTrixi.plot_bottom_triangle!(axis, reference_data.ndofs[end - 1],
                                 reference_data.ndofs[end], y_ref, polydeg + 1;
                                 trianglefontsize = 15, font = plot_font)
axislegend(axis; position = (:right, :top), labelsize = 12, font = plot_font)

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
output_path = joinpath(plots_dir, "diffusion_1d_dirichlet_penalty_convergence.pdf")
save(output_path, figure; px_per_unit = 1)
println("Saved convergence plot to: $(output_path)")
