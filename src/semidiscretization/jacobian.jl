# Select the Jacobian information supplied to SciML
abstract type AbstractJacobianStrategy end

@doc raw"""
    DefaultJacobian()

Use the default SciML Jacobian handling. HydroTrixi.jl supplies neither an analytical `jac`
function nor a `jac_prototype`.
"""
struct DefaultJacobian <: AbstractJacobianStrategy end

@doc raw"""
    SparseJacobian()

Supply a numerical sparse zero matrix as `jac_prototype`, but no analytical `jac`
function. Let ``\boldsymbol{u}_{\mathrm{state}} \in \mathbb{R}^m`` denote the state
variables supplied to the spatial operator
``\boldsymbol{\mathcal{R}}(\boldsymbol{u}_{\mathrm{state}},t)``. The prototype is composed
as

```text
spatial-operator sparsity -> physical-residual sparsity
                            -> augmented-residual sparsity
                            -> numerical residual prototype.
```

The first pattern represents
``\operatorname{pattern}(\partial\boldsymbol{\mathcal{R}}/
\partial\boldsymbol{u}_{\mathrm{state}})``. The temporal operator lifts it to the
Jacobian sparsity of the physical residual with respect to the physical state. Passive
variables ``\boldsymbol{q}``, when present, are appended last to form the augmented state
``\boldsymbol{y} =
(\boldsymbol{y}_{\mathrm{physical}},\boldsymbol{q})^\mathrm{T}`` and augmented residual
``\boldsymbol{\mathcal{F}}``. The resulting `SparseMatrixCSC` has size
`length(u_ode)` by `length(u_ode)`, has `eltype(u_ode)`, stores zeros at every structural
coordinate of
``\operatorname{pattern}(\partial\boldsymbol{\mathcal{F}}/
\partial\boldsymbol{y})``, and contains zero passive columns.

Entries arising only from the constant temporal mass matrix ``\boldsymbol{A}`` are
excluded. OrdinaryDiffEq combines the residual information with ``\boldsymbol{A}`` when
constructing the Rosenbrock matrix
``\boldsymbol{A} - \gamma\Delta t\,\boldsymbol{J}``. The algorithm's `autodiff` setting
determines how the residual Jacobian entries are computed.

The supported spatial scope is serial, nonperiodic, scalar, one-dimensional `TreeMesh`
problems using a Lobatto-Legendre `DGSEM` and `ParabolicFormulationLocalDG` with any penalty
parameter. Standard, capacity, and constitutive temporal operators and both implemented
passive-variable configurations are supported. Other signatures fail through ordinary
Julia dispatch.
"""
struct SparseJacobian <: AbstractJacobianStrategy end

# Map Jacobian strategies to SciML ODEFunction keyword arguments
function jacobian_options(::DefaultJacobian, u0_ode, semi)
    return (;)
end

function jacobian_options(::SparseJacobian, u0_ode, semi)
    jac_prototype = residual_jacobian_prototype(u0_ode, semi)
    return (; jac_prototype)
end

# Assemble the augmented residual Jacobian prototype supplied to SciML
function residual_jacobian_prototype(u_ode, semi::SemidiscretizationImplicit)
    semi_base = semi.semi_base
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi_base)
    u_state = state_variable_view(u_ode, semi)
    spatial_sparsity = spatial_operator_jacobian_sparsity(u_state, mesh, equations, solver,
                                                          semi_base.solver_parabolic, cache)
    physical_sparsity = physical_residual_jacobian_sparsity(spatial_sparsity,
                                                            semi.operator_temporal)

    return residual_jacobian_prototype(u_ode, spatial_sparsity, physical_sparsity, solver,
                                       cache, semi.passive_variables)
end

# Build the spatial-operator Jacobian sparsity for supported one-dimensional LDG schemes
function spatial_operator_jacobian_sparsity(u_state, mesh::Trixi.TreeMesh{1},
                                            equations::Trixi.AbstractEquationsParabolic{1,
                                                                                        1},
                                            solver::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                                            solver_parabolic::Trixi.ParabolicFormulationLocalDG,
                                            cache)
    # Reject MPI execution because its sparse Jacobian structure is unsupported
    if Trixi.mpi_isparallel()
        throw(ArgumentError("Sparse Jacobian structure is unsupported for MPI execution."))
    end

    # Reject periodic meshes because their sparse Jacobian structure is unsupported
    if Trixi.isperiodic(mesh)
        throw(ArgumentError("Sparse Jacobian structure is unsupported for periodic meshes."))
    end

    n_nodes = Trixi.nnodes(solver)
    n_elements = Trixi.nelements(solver, cache)
    n_interfaces = length(Trixi.eachinterface(solver, cache))
    n_state_dofs = n_nodes * n_elements
    if length(u_state) != n_state_dofs
        throw(DimensionMismatch("Spatial state has $(length(u_state)) entries, but the " *
                                "DGSEM mesh has $n_state_dofs degrees of freedom."))
    end

    expected_nonzeros = n_elements * n_nodes^2 + 2 * n_nodes * n_interfaces
    row_indices = Int[]
    column_indices = Int[]
    sizehint!(row_indices, expected_nonzeros)
    sizehint!(column_indices, expected_nonzeros)
    state_indices = LinearIndices((n_nodes, n_elements))

    # Add dense element-local spatial dependencies
    for element in Trixi.eachelement(solver, cache)
        for state_node in Trixi.eachnode(solver)
            state_dof = state_indices[state_node, element]
            for residual_node in Trixi.eachnode(solver)
                residual_dof = state_indices[residual_node, element]
                push!(row_indices, residual_dof)
                push!(column_indices, state_dof)
            end
        end
    end

    # Add the cross-element LDG interface dependencies
    for interface in Trixi.eachinterface(solver, cache)
        left_element = cache.interfaces.neighbor_ids[1, interface]
        right_element = cache.interfaces.neighbor_ids[2, interface]
        left_face_dof = state_indices[n_nodes, left_element]

        for right_node in Trixi.eachnode(solver)
            right_dof = state_indices[right_node, right_element]
            push!(row_indices, left_face_dof)
            push!(column_indices, right_dof)
            push!(row_indices, right_dof)
            push!(column_indices, left_face_dof)
        end
    end

    spatial_sparsity = sparse(row_indices, column_indices, trues(length(row_indices)),
                              n_state_dofs, n_state_dofs, |)
    if nnz(spatial_sparsity) != expected_nonzeros
        throw(ArgumentError("Spatial-operator Jacobian sparsity contains duplicate " *
                            "or inconsistent entries."))
    end

    return spatial_sparsity
end

# Retain the spatial pattern for the standard physical residual
function physical_residual_jacobian_sparsity(spatial_sparsity,
                                             ::TemporalOperatorStandard)
    return spatial_sparsity
end

# Add the nodal capacity derivative to the physical-residual sparsity
function physical_residual_jacobian_sparsity(spatial_sparsity,
                                             ::TemporalOperatorCapacity)
    n_state_dofs = size(spatial_sparsity, 1)
    identity_sparsity = sparse(1:n_state_dofs, 1:n_state_dofs, trues(n_state_dofs),
                               n_state_dofs, n_state_dofs)
    return spatial_sparsity .| identity_sparsity
end

# Embed the spatial pattern in the constitutive evolved-state block layout
function physical_residual_jacobian_sparsity(spatial_sparsity,
                                             ::TemporalOperatorConstitutive)
    n_state_dofs = size(spatial_sparsity, 1)
    zero_sparsity = spzeros(Bool, n_state_dofs, n_state_dofs)
    identity_sparsity = sparse(1:n_state_dofs, 1:n_state_dofs, trues(n_state_dofs),
                               n_state_dofs, n_state_dofs)
    return [zero_sparsity spatial_sparsity; identity_sparsity identity_sparsity]
end

# Preserve the physical-residual pattern when no passive equations are appended
function residual_jacobian_prototype(u_ode, spatial_sparsity, physical_sparsity, solver,
                                     cache, ::NoPassiveVariables)
    jac_prototype = SparseMatrixCSC{eltype(u_ode), Int}(physical_sparsity)
    fill!(nonzeros(jac_prototype), zero(eltype(jac_prototype)))
    return jac_prototype
end

# Append boundary-flux rows that depend on the adjacent state elements
function residual_jacobian_prototype(u_ode, spatial_sparsity, physical_sparsity, solver,
                                     cache, ::PassiveVariablesBoundaryFlux1D)
    n_nodes = Trixi.nnodes(solver)
    n_state_dofs = size(spatial_sparsity, 1)
    n_physical_dofs = size(physical_sparsity, 1)
    state_offset = n_physical_dofs - n_state_dofs
    state_indices = LinearIndices((n_nodes, n_state_dofs ÷ n_nodes))

    n_boundaries_per_direction = cache.boundaries.n_boundaries_per_direction
    lasts = accumulate(+, n_boundaries_per_direction)
    firsts = lasts - n_boundaries_per_direction .+ 1
    element_neg = cache.boundaries.neighbor_ids[firsts[1]]
    element_pos = cache.boundaries.neighbor_ids[firsts[2]]

    boundary_sparsity = spzeros(Bool, 2, n_physical_dofs)
    for state_node in Trixi.eachnode(solver)
        boundary_sparsity[1, state_offset + state_indices[state_node, element_neg]] = true
        boundary_sparsity[2, state_offset + state_indices[state_node, element_pos]] = true
    end

    passive_columns = spzeros(Bool, n_physical_dofs + 2, 2)
    augmented_sparsity = [physical_sparsity; boundary_sparsity]
    augmented_sparsity = [augmented_sparsity passive_columns]
    jac_prototype = SparseMatrixCSC{eltype(u_ode), Int}(augmented_sparsity)
    fill!(nonzeros(jac_prototype), zero(eltype(jac_prototype)))
    return jac_prototype
end
