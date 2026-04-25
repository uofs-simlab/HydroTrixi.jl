using DiffuSEM
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

function run_level(level)
    case = mixed_dirichlet_neumann_case(;
                                        diffusivity = diffusivity,
                                        forcing_amplitude = forcing_amplitude,
                                        forcing_frequency = forcing_frequency,
                                        dirichlet_mean = dirichlet_mean,
                                        tspan = tspan,
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

ndofs = Int[]
l2_errors = Float64[]
linf_errors = Float64[]
dts = Float64[]

for level in levels
    dofs, l2, linf, dt = run_level(level)
    push!(ndofs, dofs)
    push!(l2_errors, l2)
    push!(linf_errors, linf)
    push!(dts, dt)
end

l2_eoc = compute_eoc(l2_errors)
linf_eoc = compute_eoc(linf_errors)

println("1D diffusion equation convergence study (mixed Dirichlet-Neumann, LDG)")
println("diffusivity = $(diffusivity), polydeg = $(polydeg), tspan = $(tspan)")
println("forcing_amplitude = $(forcing_amplitude), forcing_frequency = $(forcing_frequency)")
println()
@printf("%-7s %-10s %-11s %-14s %-8s %-14s %-8s\n",
        "level",
        "ndofs",
        "dt",
        "L2 error",
        "EOC",
        "Linf error",
        "EOC")
for i in eachindex(levels)
    level = levels[i]
    @printf("%-7d %-10d %-11.3e %-14.6e ", level, ndofs[i], dts[i], l2_errors[i])
    isnan(l2_eoc[i]) ? @printf("%-8s ", "-") : @printf("%-8.3f ", l2_eoc[i])
    @printf("%-14.6e ", linf_errors[i])
    isnan(linf_eoc[i]) ? @printf("%-8s\n", "-") : @printf("%-8.3f\n", linf_eoc[i])
end

plot_convergence_1d(ndofs,
                    l2_errors,
                    linf_errors;
                    output_path = joinpath(plots_dir,
                                           "diffusion_equation_1d_mixed_dirichlet_neumann_convergence.pdf"),
                    triangle_order = polydeg + 1,)
