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

schemes = ((; name = "LDG", penalty_factor = 0,
            labels = (L"\mathrm{LDG}\ L^2", L"\mathrm{LDG}\ L^\infty"), color = 1),
           (; name = "LDG + boundary penalty", penalty_factor = 1,
            labels = (L"\mathrm{LDG + boundary\ penalty}\ L^2",
                      L"\mathrm{LDG + boundary\ penalty}\ L^\infty"),
            color = 2))

function run_scheme(scheme)
    _, errors = convergence_test(@__MODULE__, elixir, iterations;
                                 diffusivity = diffusivity, polydeg = polydeg,
                                 tspan = tspan,
                                 initial_refinement_level = base_initial_refinement_level,
                                 solver_parabolic = ParabolicFormulationLocalDG(),
                                 penalty_factor = scheme.penalty_factor)

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
series = map(schemes) do scheme
    data = results[scheme.name]
    (; x = data.ndofs, errors = (data.l2_errors, data.linf_errors),
     labels = scheme.labels, color = scheme.color)
end

plots_dir = mkpath(joinpath(dirname(dirname(@__DIR__)), "plots"))
output_path = joinpath(plots_dir, "diffusion_1d_dirichlet_penalty_convergence.pdf")
plot_convergence_1d(series; output_path, triangle_order = polydeg + 1,
                    trianglefontsize = 15, legendfontsize = 12)
println("Saved convergence plot to: $(output_path)")
