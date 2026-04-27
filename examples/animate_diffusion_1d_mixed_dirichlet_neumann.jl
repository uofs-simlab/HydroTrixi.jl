# Example: save a time-resolved animation for an existing mixed
# Dirichlet-Neumann diffusion solution.

using DiffuSEM
plots_dir = mkpath(joinpath(dirname(@__DIR__), "plots"))
if !@isdefined(animation_format)
    animation_format = "mp4"
end

output_path = joinpath(plots_dir,
                       "diffusion_equation_1d_mixed_dirichlet_neumann_solution_long.$(animation_format)")

animate_solution_1d(sol;
                    exact_solution = exact_solution,
                    output_path = output_path,
                    ylims = (0.5, 1.5),
                    framerate = 30)

println("Saved solution animation to: $(output_path)")
@show sol.t[end]
