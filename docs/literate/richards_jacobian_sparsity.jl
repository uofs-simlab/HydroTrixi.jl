# ```@meta
# CurrentModule = HydroTrixi
# ```
#
# # Sparse Jacobian evaluation for implicit semidiscretizations
#
# `SemidiscretizationImplicit` supplies a sparse residual Jacobian prototype to
# OrdinaryDiffEq.jl by default. The prototype determines the Jacobian storage, while the
# time integration algorithm's `autodiff` setting determines how its entries are computed.
# With adaptive mesh refinement, HydroTrixi.jl rebuilds the prototype and associated
# Rosenbrock solver state whenever the mesh topology changes.
#
# For the mixed Richards formulation, the state is
# $\boldsymbol{y}=(\boldsymbol{\Theta},\boldsymbol{\Psi})^\mathrm{T}$, where
# $\boldsymbol{\Theta}$ is water content, and $\boldsymbol{\Psi}$ is pressure head. Let
# $\boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t)$ denote the spatial residual, and let
# $\boldsymbol{\vartheta}(\boldsymbol{\Psi})$ denote the nodal constitutive map. The
# residual and its Jacobian sparsity pattern are
#
# ```math
# \boldsymbol{\mathcal{F}}(\boldsymbol{y},t) =
# \begin{bmatrix}
# \boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t) \\
# \boldsymbol{\Theta}-\boldsymbol{\vartheta}(\boldsymbol{\Psi})
# \end{bmatrix},
# \qquad
# \operatorname{pattern}\left(\frac{\partial\boldsymbol{\mathcal{F}}}
# {\partial\boldsymbol{y}}\right)
# =
# \begin{bmatrix}
# \boldsymbol{0} & \boldsymbol{S} \\
# \boldsymbol{I} & \boldsymbol{I}
# \end{bmatrix},
# \qquad
# \boldsymbol{S} = \operatorname{pattern}\left(
# \frac{\partial\boldsymbol{\mathcal{R}}}{\partial\boldsymbol{\Psi}}\right).
# ```
#
# `jac_prototype` stores zeros at the coordinates in this pattern. Entries arising only
# from the temporal mass matrix $\boldsymbol{A}$ are added by OrdinaryDiffEq.jl when it
# constructs the Rosenbrock matrix $\boldsymbol{A}-\gamma\Delta t_n\boldsymbol{J}_n$.
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

# For this one-dimensional mixed Richards discretization, the residual prototype has
# $K(N+1)^2 + 2(N+1)(K-1) + 2K(N+1)$ stored entries. The first term contains dense
# element-local entries of
# $\partial_{\boldsymbol{\Psi}}\boldsymbol{\mathcal{R}}$, the second term contains
# LDG interface entries, and the final term contains the two constitutive identity
# diagonals. We verify the prototype size, number of stored entries, and stored values:

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

# Next, we visualize the residual Jacobian sparsity pattern:

interface_pattern = spzeros(Bool, 2 * n_state_dofs, 2 * n_state_dofs)
node_indices = LinearIndices((n_nodes, n_elements))
for interface in Trixi.eachinterface(solver, cache)
    left_element = cache.interfaces.neighbor_ids[1, interface]
    right_element = cache.interfaces.neighbor_ids[2, interface]
    left_face_dof = node_indices[n_nodes, left_element]
    for right_node in Trixi.eachnode(solver)
        right_dof = node_indices[right_node, right_element]
        interface_pattern[left_face_dof, n_state_dofs + right_dof] = true
        interface_pattern[right_dof, n_state_dofs + left_face_dof] = true
    end
end

figure = Figure(size = (1050, 500))
jacobian_axis = Axis(figure[1, 1]; aspect = DataAspect(), yreversed = true)
spatial_axis = Axis(figure[1, 2]; aspect = DataAspect(), yreversed = true)
hidedecorations!(jacobian_axis)
hidedecorations!(spatial_axis)
spy!(jacobian_axis, transpose(jac_prototype); color = :black)
spy!(jacobian_axis, transpose(interface_pattern); color = :orange)
spy!(spatial_axis,
     transpose(jac_prototype[1:n_state_dofs,
                             (n_state_dofs + 1):(2 * n_state_dofs)]); color = :black)
spy!(spatial_axis,
     transpose(interface_pattern[1:n_state_dofs,
                                 (n_state_dofs + 1):(2 * n_state_dofs)]);
     color = :orange)
vlines!(jacobian_axis, 1:(2 * n_state_dofs - 1); color = :white, linewidth = 0.5)
hlines!(jacobian_axis, 1:(2 * n_state_dofs - 1); color = :white, linewidth = 0.5)
vlines!(spatial_axis, 1:(n_state_dofs - 1); color = :white, linewidth = 0.5)
hlines!(spatial_axis, 1:(n_state_dofs - 1); color = :white, linewidth = 0.5)
vlines!(jacobian_axis, n_state_dofs; color = :black, linewidth = 1)
hlines!(jacobian_axis, n_state_dofs; color = :black, linewidth = 1)

figure_path = joinpath(asset_dir, "pattern.png")
save(figure_path, figure; px_per_unit = 2)
println("Saved sparse Jacobian pattern with $(nnz(jac_prototype)) entries to " *
        figure_path)

# In the figure below, the residual Jacobian sparsity is shown on the left, and the
# spatial-operator Jacobian sparsity
# $\partial_{\boldsymbol{\Psi}}\boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t)$ is shown on
# the right. The orange entries show dependencies across element interfaces introduced by
# the LDG fluxes.
#
# ![Sparse Jacobian and spatial block patterns](../assets/generated/richards_jacobian_sparsity/pattern.png)
