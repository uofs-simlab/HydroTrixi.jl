using CairoMakie
using DiffuSEM
using LaTeXStrings
using Printf: @printf
using Trixi

include(joinpath(dirname(@__FILE__), "diffusion_equation_1d_mixed_dirichlet_neumann.jl"))

plots_dir = joinpath(dirname(@__DIR__), "plots")
mkpath(plots_dir)

diffusivity = 0.5
forcing_amplitude = 0.4
forcing_frequency = 4.0
dirichlet_mean = 1.0
polydeg = 3
tspan = (0.0, 0.25)
levels = collect(2:5)
dt_factor = 0.01

schemes = ((; name = "LDG", solver_parabolic = ParabolicFormulationLocalDG()),
           (; name = "BR1", solver_parabolic = ParabolicFormulationBassiRebay1()))

function run_level(level, scheme)
    case = mixed_dirichlet_neumann_case(;
                                        diffusivity = diffusivity,
                                        forcing_amplitude = forcing_amplitude,
                                        forcing_frequency = forcing_frequency,
                                        dirichlet_mean = dirichlet_mean,
                                        tspan = tspan,
                                        solver_parabolic = scheme.solver_parabolic,
                                        initial_refinement_level = level,
                                        polydeg = polydeg,)

    h = 1.0 / (2.0^level)
    dt = dt_factor * h^2
    sol = solve_mixed_dirichlet_neumann_case(case;
                                             dt = dt,
                                             adaptive = false,
                                             save_everystep = false,)

    analysis_callback = AnalysisCallback(case.semi, interval = typemax(Int))
    errors = analysis_callback(sol)

    return length(sol.u[end]), errors.l2[1], errors.linf[1], dt
end

function run_scheme(scheme)
    ndofs = Int[]
    l2_errors = Float64[]
    linf_errors = Float64[]
    dts = Float64[]

    for level in levels
        dofs, l2, linf, dt = run_level(level, scheme)
        push!(ndofs, dofs)
        push!(l2_errors, l2)
        push!(linf_errors, linf)
        push!(dts, dt)
    end

    return (;
            ndofs,
            l2_errors,
            linf_errors,
            dts,
            l2_eoc = compute_eoc(l2_errors),
            linf_eoc = compute_eoc(linf_errors),)
end

results = Dict(scheme.name => run_scheme(scheme) for scheme in schemes)

println("1D diffusion equation convergence study (mixed Dirichlet-Neumann)")
println("Comparison: BR1 vs LDG parabolic schemes")
println("diffusivity = $(diffusivity), polydeg = $(polydeg), tspan = $(tspan)")
println("forcing_amplitude = $(forcing_amplitude), forcing_frequency = $(forcing_frequency)")
println()

for scheme in schemes
    data = results[scheme.name]
    println("Scheme: $(scheme.name)")
    @printf("%-7s %-10s %-11s %-14s %-8s %-14s %-8s\n",
            "level",
            "ndofs",
            "dt",
            "L2 error",
            "EOC",
            "Linf error",
            "EOC")
    for i in eachindex(levels)
        @printf("%-7d %-10d %-11.3e %-14.6e ",
                levels[i],
                data.ndofs[i],
                data.dts[i],
                data.l2_errors[i])
        isnan(data.l2_eoc[i]) ? @printf("%-8s ", "-") : @printf("%-8.3f ", data.l2_eoc[i])
        @printf("%-14.6e ", data.linf_errors[i])
        isnan(data.linf_eoc[i]) ? @printf("%-8s\n", "-") :
        @printf("%-8.3f\n", data.linf_eoc[i])
    end
    println()
end

ldg = results["LDG"]
br1 = results["BR1"]

println("BR1/LDG error ratio by level")
@printf("%-7s %-12s %-12s\n", "level", "L2 ratio", "Linf ratio")
for i in eachindex(levels)
    @printf("%-7d %-12.4f %-12.4f\n",
            levels[i],
            br1.l2_errors[i]/ldg.l2_errors[i],
            br1.linf_errors[i]/ldg.linf_errors[i])
end

DiffuSEM.set_serif_tex_theme!()
plot_font = DiffuSEM.DEFAULT_PLOT_FONT

fig = Figure(size = DiffuSEM.DEFAULT_CONVERGENCE_FIGSIZE, fontsize = 15)
xticks = DiffuSEM._doubling_dof_ticks(ldg.ndofs; base = 10)
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
              linewidth = 1.8,)
scatterlines!(ax,
              br1.ndofs,
              br1.l2_errors;
              label = L"\mathrm{BR1}\ L^2",
              color = colors[2],
              linestyle = :solid,
              linewidth = 1.8,)
scatterlines!(ax,
              ldg.ndofs,
              ldg.linf_errors;
              label = L"\mathrm{LDG}\ L^\infty",
              color = colors[1],
              linestyle = :dash,
              linewidth = 1.8,)
scatterlines!(ax,
              br1.ndofs,
              br1.linf_errors;
              label = L"\mathrm{BR1}\ L^\infty",
              color = colors[2],
              linestyle = :dash,
              linewidth = 1.8,)

axislegend(ax; position = (:right, :top), labelsize = 14, font = plot_font)

comparison_plot = joinpath(plots_dir,
                           "diffusion_equation_1d_mixed_dirichlet_neumann_convergence_br1_vs_ldg.pdf")
save(comparison_plot, fig; px_per_unit = 1)
comparison_plot_png = joinpath(plots_dir,
                               "diffusion_equation_1d_mixed_dirichlet_neumann_convergence_br1_vs_ldg.png")
save(comparison_plot_png, fig; px_per_unit = 1)
println()
println("Saved comparison plot to: $(comparison_plot)")
println("Saved comparison plot to: $(comparison_plot_png)")
