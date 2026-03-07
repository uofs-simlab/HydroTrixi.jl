using DiffuSEM
include(joinpath(dirname(@__FILE__), "diffusion_equation_1d_mixed_dirichlet_neumann.jl"))

plots_dir = joinpath(dirname(@__DIR__), "plots")
mkpath(plots_dir)

case = mixed_dirichlet_neumann_case()
sol = solve_mixed_dirichlet_neumann_case(case;
                                         dt=5.0e-5,
                                         adaptive=false,
                                         save_everystep=false)

plot_solution_1d(sol;
                 exact_solution=case.exact_solution,
                 output_path=joinpath(plots_dir,
                                      "diffusion_equation_1d_mixed_dirichlet_neumann_solution.pdf"))

@show sol.t[end]
