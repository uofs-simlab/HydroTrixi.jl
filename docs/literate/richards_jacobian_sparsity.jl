# ```@meta
# CurrentModule = HydroTrixi
# ```
#
# # Sparse Jacobian evaluation for the mixed-form Richards equation
#
# HydroTrixi.jl can supply a known Jacobian sparsity pattern to OrdinaryDiffEq.jl so that
# Jacobian evaluation and the associated linear algebraic operations can exploit sparse
# matrices. The time-integration algorithm's `autodiff` setting determines how Jacobian
# entries are computed. The default `AutoFiniteDiff()` setting in HydroTrixi.jl uses
# graph-coloured finite differences with the sparse prototype. This is useful for
# implicit time integration, particularly for Rosenbrock-Wanner methods that construct
# linear systems from an exact or approximate Jacobian matrix. You can enable this path as
# follows when creating the ODE problem:
#
# ```julia
# ode = semidiscretize(semi, tspan; jacobian = SparseJacobian())
# ```
#
# This option supplies a structural prototype for Jacobian construction and linear algebra.
# The default `jacobian = DefaultJacobian()` supplies no Jacobian information, and 
# OrdinaryDiffEq.jl must resort to the most general path, which treats the Jacobian as a 
# dense matrix and therefore requires far more function evaluations and memory than 
# necessary for a sparse Jacobian such as that resulting from an LDG discretization of the 
# Richards equation. With adaptive mesh refinement, HydroTrixi.jl rebuilds the sparse 
# prototype and the associated Rosenbrock solver state whenever the mesh topology changes.
#
# For $K$ mesh elements and polynomial degree $N$, let $\boldsymbol{\Theta}(t)$ and
# $\boldsymbol{\Psi}(t)$ denote the global water-content and pressure-head vectors,
# respectively. Each vector has $K(N+1)$ entries. Consistent with the mixed formulation,
# the global state, mass matrix, and nonlinear right-hand side are
#
# ```math
# \boldsymbol{y}(t) =
# \begin{bmatrix}
# \boldsymbol{\Theta}(t) \\ \boldsymbol{\Psi}(t)
# \end{bmatrix},
# \qquad
# \boldsymbol{A} =
# \begin{bmatrix}
# \boldsymbol{I} & \boldsymbol{0} \\
# \boldsymbol{0} & \boldsymbol{0}
# \end{bmatrix},
# \qquad
# \boldsymbol{\mathcal{F}}(\boldsymbol{y}(t),t) =
# \begin{bmatrix}
# \boldsymbol{\mathcal{R}}(\boldsymbol{\Psi}(t),t) \\
# \boldsymbol{\Theta}(t) -
# \boldsymbol{\vartheta}(\boldsymbol{\Psi}(t))
# \end{bmatrix},
# ```
#
# where $\boldsymbol{\mathcal{R}}$ is the global LDG spatial residual and
# $\boldsymbol{\vartheta}$ applies the water-content function at each node. This is the
# mass-matrix system
# $\boldsymbol{A}\dot{\boldsymbol{y}} =
# \boldsymbol{\mathcal{F}}(\boldsymbol{y},t)$ from the mixed formulation. The residual
# Jacobian has the block structure
#
# ```math
# \boldsymbol{J}(\boldsymbol{y},t)
# = \frac{\partial\boldsymbol{\mathcal{F}}}{\partial\boldsymbol{y}}(\boldsymbol{y},t)
# =
# \begin{bmatrix}
# \boldsymbol{0} &
# \partial_{\boldsymbol{\Psi}}
# \boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t) \\
# \boldsymbol{I} & -\boldsymbol{C}(\boldsymbol{\Psi})
# \end{bmatrix}.
# ```
#
# Here, $\boldsymbol{C}(\boldsymbol{\Psi})$ is the diagonal matrix of nodal capacity
# values. HydroTrixi.jl supplies the corresponding sparse residual pattern as the
# `jac_prototype` used by OrdinaryDiffEq.jl. For a Rosenbrock-Wanner step, OrdinaryDiffEq
# solves systems with $\boldsymbol{A} - \gamma \Delta t_n \boldsymbol{J}_n$ and augments
# its solver matrices with the differential-block diagonal from $\boldsymbol{A}$.
#
# We now construct a small problem with $K=4$ and $N=3$ following
# `examples/elixir_richards_manufactured_solution.jl`.

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
mesh = TreeMesh(problem.domain...; initial_refinement_level = 2, n_cells_max = 100,
                periodicity = false)
solver = DGSEM(polydeg = 3, surface_flux = flux_central)
semi = SemidiscretizationImplicit(mesh, problem, solver;
                                  solver_parabolic = ParabolicFormulationLocalDG(),
                                  form = MixedForm(),
                                  passive_variables = NoPassiveVariables())

ode = semidiscretize(semi, problem.tspan; jacobian = SparseJacobian());

# For this one-dimensional Richards discretization, the Jacobian has
# $K(N+1)^2 + 2(N+1)(K-1) + 2K(N+1)$ stored entries. The first term contains dense
# element-local entries of
# $\partial_{\boldsymbol{\Psi}}\boldsymbol{\mathcal{R}}$, the second term contains
# LDG interface entries, and the final term contains the diagonal capacity entries of
# $-\boldsymbol{C}(\boldsymbol{\Psi})$ and the identity entries of
# $\partial_{\boldsymbol{\Theta}}\boldsymbol{\mathcal{F}}$. We verify the sparsity
# pattern size and number of nonzeros against this formula:

mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
n_nodes = Trixi.nnodes(solver)
n_elements = Trixi.nelements(solver, cache)
n_interfaces = length(Trixi.eachinterface(solver, cache))
n_physical_dofs = n_nodes * n_elements
expected_nonzeros = n_elements * n_nodes^2 + 2 * n_nodes * n_interfaces +
                    2 * n_physical_dofs

@assert size(ode.f.sparsity) == (2 * n_physical_dofs, 2 * n_physical_dofs)
@assert nnz(ode.f.sparsity) == expected_nonzeros

# Next, we will visualize the sparsity pattern of the Jacobian.

interface_sparsity = spzeros(Bool, 2 * n_physical_dofs, 2 * n_physical_dofs)
physical_indices = LinearIndices((n_nodes, n_elements))
for interface in Trixi.eachinterface(solver, cache)
    left_element = cache.interfaces.neighbor_ids[1, interface]
    right_element = cache.interfaces.neighbor_ids[2, interface]
    left_face_dof = physical_indices[n_nodes, left_element]
    for right_node in Trixi.eachnode(solver)
        right_dof = physical_indices[right_node, right_element]
        interface_sparsity[left_face_dof, n_physical_dofs + right_dof] = true
        interface_sparsity[right_dof, n_physical_dofs + left_face_dof] = true
    end
end

figure = Figure(size = (1050, 500))
jacobian_axis = Axis(figure[1, 1]; aspect = DataAspect(), yreversed = true)
spatial_axis = Axis(figure[1, 2]; aspect = DataAspect(), yreversed = true)
hidedecorations!(jacobian_axis)
hidedecorations!(spatial_axis)
spy!(jacobian_axis, transpose(ode.f.sparsity); color = :black)
spy!(jacobian_axis, transpose(interface_sparsity); color = :orange)
spy!(spatial_axis,
     transpose(ode.f.sparsity[1:n_physical_dofs,
                              (n_physical_dofs + 1):(2 * n_physical_dofs)]); color = :black)
spy!(spatial_axis,
     transpose(interface_sparsity[1:n_physical_dofs,
                                  (n_physical_dofs + 1):(2 * n_physical_dofs)]);
     color = :orange)
vlines!(jacobian_axis, 1:(2 * n_physical_dofs - 1); color = :white, linewidth = 0.5)
hlines!(jacobian_axis, 1:(2 * n_physical_dofs - 1); color = :white, linewidth = 0.5)
vlines!(spatial_axis, 1:(n_physical_dofs - 1); color = :white, linewidth = 0.5)
hlines!(spatial_axis, 1:(n_physical_dofs - 1); color = :white, linewidth = 0.5)
vlines!(jacobian_axis, n_physical_dofs; color = :black, linewidth = 1)
hlines!(jacobian_axis, n_physical_dofs; color = :black, linewidth = 1)

figure_path = joinpath(asset_dir, "pattern.png")
save(figure_path, figure; px_per_unit = 2)
println("Saved sparse Jacobian pattern with $(nnz(ode.f.sparsity)) entries to " *
        figure_path)

# In the figure below, the sparsity pattern of the full residual Jacobian
# $\boldsymbol{J}(\boldsymbol{y},t)$ is shown on the left, and the sparsity pattern of the
# spatial block
# $\partial_{\boldsymbol{\Psi}}\boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t)$ is shown on
# the right. The orange entries show dependencies across element interfaces introduced by
# the LDG fluxes.
#
# ![Sparse Jacobian and spatial block patterns](../assets/generated/richards_jacobian_sparsity/pattern.png)
