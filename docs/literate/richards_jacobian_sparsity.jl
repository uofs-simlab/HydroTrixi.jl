# ```@meta
# CurrentModule = HydroTrixi
# ```
#
# # Sparse Jacobian evaluation for implicit semidiscretizations
#
# By default, `semidiscretize` supplies a sparse Jacobian prototype and a column
# colouring to OrdinaryDiffEq.jl. The prototype specifies which matrix entries are stored.
# The time integration algorithm's `autodiff` option selects how their values are computed.
#
# In mixed form, the state is
# $\boldsymbol{y}=(\boldsymbol{\Theta},\boldsymbol{\Psi})^\mathrm{T}$, where
# $\boldsymbol{\Theta}$ contains nodal water contents and $\boldsymbol{\Psi}$ contains
# pressure heads. Let $\boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t)$ denote the spatial
# residual and $\boldsymbol{\vartheta}(\boldsymbol{\Psi})$ the water contents evaluated
# from pressure head using the constitutive model. The residual and its Jacobian sparsity
# pattern are
#
# ```math
# \boldsymbol{\mathcal{F}}(\boldsymbol{y},t) =
# \begin{bmatrix}
# \boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t) \\
# \boldsymbol{\Theta}-\boldsymbol{\vartheta}(\boldsymbol{\Psi})
# \end{bmatrix},
# \quad
# \operatorname{pattern}\left(\frac{\partial\boldsymbol{\mathcal{F}}}
# {\partial\boldsymbol{y}}\right)
# =
# \begin{bmatrix}
# \boldsymbol{0} & \boldsymbol{S} \\
# \boldsymbol{I} & \boldsymbol{I}
# \end{bmatrix},
# \quad
# \boldsymbol{S} = \operatorname{pattern}\left(
# \frac{\partial\boldsymbol{\mathcal{R}}}{\partial\boldsymbol{\Psi}}\right).
# ```
#
# The lower-right $\boldsymbol{I}$ denotes a diagonal sparsity pattern. The corresponding
# Jacobian block contains the negative nodal water capacities.
#
# `jac_prototype` stores explicit zeros at the positions allowed by this pattern.
# OrdinaryDiffEq.jl fills these positions with derivative values and combines the Jacobian
# with the mass matrix $\boldsymbol{M}$ to form the Rosenbrock matrix:
#
# ```math
# \boldsymbol{M}-\gamma\Delta t^n
# \frac{\partial\boldsymbol{\mathcal{F}}}{\partial\boldsymbol{y}}
# (\boldsymbol{y}^n,t^n).
# ```
#
# We now construct and visualize this pattern for a small mixed Richards problem with $K=4$
# elements and polynomial degree $N=3$, following
# `examples/elixirs/elixir_richards_manufactured_solution.jl`:

using CairoMakie
using HydroTrixi
using SparseArrays
using Trixi

tutorial_utils_root = get(ENV, "HYDROTRIXI_DOCS_LITERATE", @__DIR__) #hide
tutorial_utils_path = joinpath(tutorial_utils_root, "tutorial_utils.jl") #hide
include(tutorial_utils_path) #hide
using .TutorialUtils: docs_generated_dir #hide
nothing #hide

asset_dir = docs_generated_dir("richards_jacobian_sparsity")

problem = HydrologicProblemRichardsManufacturedSolution(tspan = (0.0, 1.0))
mesh = TreeMesh(problem.domain...; initial_refinement_level = 2, periodicity = false)
solver = DGSEM(polydeg = 3)
semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  form = MixedForm(),
                                  passive_variables = NoPassiveVariables())

ode = semidiscretize(semi, problem.tspan);

# For this one-dimensional mixed Richards discretization, the Jacobian prototype has
# $K(N+1)^2 + 2(N+1)(K-1) + 2K(N+1)$ stored entries. The three terms count the
# element-local spatial couplings, the couplings across element interfaces, and the two
# diagonal blocks from the constitutive constraint. We check the matrix dimensions,
# stored-entry count, and zero initialization below.

mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
n_nodes = Trixi.nnodes(solver)
n_elements = Trixi.nelements(solver, cache)
n_interfaces = length(Trixi.eachinterface(solver, cache))
n_state_dofs = n_nodes * n_elements
expected_nonzeros = n_elements * n_nodes^2 + 2 * n_nodes * n_interfaces +
                    2 * n_state_dofs
jac_prototype = ode.f.jac_prototype

@assert size(jac_prototype) == (2 * n_state_dofs, 2 * n_state_dofs)
@assert nnz(jac_prototype) == expected_nonzeros
@assert all(iszero, nonzeros(jac_prototype))

# ## Column colouring
#
# Columns can share a colour when they have no stored entries in the same row. Their
# derivatives can then be recovered from a single differentiation direction. HydroTrixi.jl
# supplies the colour assigned to each column in `ode.f.colorvec`.

column_colors = ode.f.colorvec
n_spatial_colors = 2 * n_nodes + 1
n_colors = maximum(column_colors)

@assert column_colors[1:n_state_dofs] == fill(n_colors, n_state_dofs)
@assert column_colors[(n_state_dofs + 1):end] ==
        mod1.(1:n_state_dofs, n_spatial_colors)
@assert n_colors == n_spatial_colors + 1

# The pressure-head columns use $2N+3$ colours, repeated in storage order. The water-content
# columns share one additional colour because their entries occupy distinct rows. We check
# that each row contains at most one stored entry of each colour.

row_indices, column_indices, _ = findnz(jac_prototype)
for row in axes(jac_prototype, 1)
    row_colors = column_colors[column_indices[row_indices .== row]]
    @assert allunique(row_colors)
end

# To visualize the sparsity pattern and colouring together, we assign each stored Jacobian
# entry the colour of its column. The black lines separate the water-content and
# pressure-head blocks.

colored_pattern = fill(NaN, size(jac_prototype))
for (row, column) in zip(row_indices, column_indices)
    colored_pattern[row, column] = column_colors[column]
end

palette = cgrad(:tab10, n_colors; categorical = true)
figure = Figure(size = (650, 650))
jacobian_axis = Axis(figure[1, 1]; aspect = DataAspect(), yreversed = true)
hidedecorations!(jacobian_axis)
heatmap!(jacobian_axis, transpose(colored_pattern); colormap = palette,
         colorrange = (0.5, n_colors + 0.5))
vlines!(jacobian_axis, n_state_dofs + 0.5; color = :black, linewidth = 1)
hlines!(jacobian_axis, n_state_dofs + 0.5; color = :black, linewidth = 1)

figure_path = joinpath(asset_dir, "coloring.png")
save(figure_path, figure; px_per_unit = 2)
println("Saved $(n_colors)-colour Jacobian colouring with $(nnz(jac_prototype)) entries " *
        "to " * figure_path)

# ![Sparse Jacobian column colouring](../assets/generated/richards_jacobian_sparsity/coloring.png)
#
# !!! note "Adaptive mesh refinement"
#     HydroTrixi.jl rebuilds the prototype and associated Rosenbrock solver state whenever
#     the mesh topology changes.
