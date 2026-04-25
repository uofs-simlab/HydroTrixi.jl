using DiffuSEM
include(joinpath(dirname(@__FILE__), "diffusion_equation_1d_mixed_dirichlet_neumann.jl"))

plots_dir = joinpath(dirname(@__DIR__), "plots")
mkpath(plots_dir)

case = mixed_dirichlet_neumann_case()
sol = solve_mixed_dirichlet_neumann_case(case;
                                         dt = 5.0e-5,
                                         adaptive = false,
                                         saveat = range(case.tspan[1], case.tspan[2],
                                                        length = 301),)

animate_solution_1d(sol;
                    exact_solution = case.exact_solution,
                    output_path = joinpath(plots_dir,
                                           "diffusion_equation_1d_mixed_dirichlet_neumann_solution_long.mp4"),
                    ylims = (0.55, 1.45),
                    framerate = 30,)

@show sol.t[end]
