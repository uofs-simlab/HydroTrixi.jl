# Select the Jacobian information supplied to SciML
@doc raw"""
    AbstractJacobianStrategy

Abstract supertype for strategies that control the Jacobian information supplied by
[`semidiscretize`](@ref) for a [`SemidiscretizationImplicit`](@ref). The Jacobian strategy
is specified by the `jacobian` keyword argument to [`semidiscretize`](@ref) and defines
whether a sparse Jacobian prototype is supplied (the default `SparseJacobian()` option) or
whether a dense Jacobian is used (`DenseJacobian()`). This is independent of the
differentiation backend, which controls how the Jacobian entries are computed and is
specified by the time integration algorithm's `autodiff` keyword argument.
"""
abstract type AbstractJacobianStrategy end

@doc raw"""
    DenseJacobian()

When [`semidiscretize`](@ref) is called with `jacobian = DenseJacobian()`, HydroTrixi.jl
selects dense Jacobian storage and does not supply an analytical `jac` function.

The Jacobian entries are computed using the backend passed to the
time integration algorithm's `autodiff` keyword, for example, `ADTypes.AutoForwardDiff()`
or `ADTypes.AutoFiniteDiff()`.
"""
struct DenseJacobian <: AbstractJacobianStrategy end

@doc raw"""
    SparseJacobian()

When [`semidiscretize`](@ref) is called with `jacobian = SparseJacobian()`, HydroTrixi.jl
supplies a sparse zero matrix as `jac_prototype`, but no analytical `jac` function. The
resulting `SparseMatrixCSC` is a prototype that encodes the full Jacobian sparsity pattern.
It has size `length(u_ode)` by `length(u_ode)`, has `eltype(u_ode)`, and stores zeros at
every coordinate in the pattern. When passive variables are present, their rows are
included, and their columns are zero.

Entries arising only from the constant temporal mass matrix ``\boldsymbol{A}`` are
excluded, since OrdinaryDiffEq.jl combines the Jacobian information with
``\boldsymbol{A}`` when constructing the Rosenbrock matrix
``\boldsymbol{A} - \gamma\Delta t\,\boldsymbol{J}``.

The Jacobian entries are computed using the backend passed to the
time integration algorithm's `autodiff` keyword, for example, `ADTypes.AutoForwardDiff()`
or `ADTypes.AutoFiniteDiff()`.

This is the default Jacobian strategy for `SemidiscretizationImplicit`. Jacobian storage,
the time integration algorithm, and the differentiation backend are configured separately.

For the supported serial, nonperiodic, one-dimensional LGL-DGSEM LDG discretization,
HydroTrixi.jl also supplies a deterministic analytical column colouring. The construction
assumes a fixed polynomial degree during h-adaptive mesh refinement and Trixi.jl's current
physical ordering of one-dimensional leaf cells. It must be revisited if the stencil,
field coupling, mesh type, boundary treatment, or leaf ordering changes.
"""
struct SparseJacobian <: AbstractJacobianStrategy end

# Map Jacobian strategies to SciML ODEFunction keyword arguments
function jacobian_options(::DenseJacobian, u0_ode, semi)
    return (;)
end

function jacobian_options(::SparseJacobian, u0_ode, semi)
    jac_prototype = jacobian_prototype(u0_ode, semi)
    return (; jac_prototype, colorvec = colorvec(u0_ode, semi))
end

# Construct a deterministic column colouring for the supported one-dimensional LDG
# pattern. With N + 1 nodes per element, every largest spatial conflict set contains the
# right endpoint of one element and every node of the next two elements, requiring
# 2N + 3 = 2 * n_nodes + 1 colours. The contiguous Trixi.jl element layout lets us assign
# these colours cyclically in storage order.
function colorvec(u_ode, semi::SemidiscretizationImplicit)
    n_spatial_colors = 2 * Trixi.nnodes(semi.semi_base.solver) + 1

    # SparseMatrixColorings requires every colour label to be no larger than the number of
    # columns. This cap matters only for exceptionally small systems that cannot contain
    # the full conflict set.
    auxiliary_color = min(n_spatial_colors + 1, length(u_ode))
    colors = fill(auxiliary_color, length(u_ode))
    spatial_colors = state_variable_view(colors, semi)
    spatial_colors .= mod1.(eachindex(spatial_colors), n_spatial_colors)
    return colors
end

# Assemble the Jacobian prototype supplied to SciML
function jacobian_prototype(u_ode, semi::SemidiscretizationImplicit)
    semi_base = semi.semi_base
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi_base)
    u_state = state_variable_view(u_ode, semi)
    spatial_pattern = spatial_operator_jacobian_sparsity_pattern(u_state, mesh, equations,
                                                                 solver,
                                                                 semi_base.solver_parabolic,
                                                                 cache)
    physical_pattern = physical_jacobian_sparsity_pattern(spatial_pattern,
                                                          semi.operator_temporal)

    return jacobian_prototype(u_ode, spatial_pattern, physical_pattern, solver, cache,
                              semi.passive_variables)
end

# Build the spatial operator Jacobian pattern for supported one-dimensional LDG schemes
function spatial_operator_jacobian_sparsity_pattern(u_state, mesh::Trixi.TreeMesh{1},
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
        throw(ArgumentError("Sparse Jacobian structure is unsupported for " *
                            "periodic meshes."))
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

    spatial_pattern = sparse(row_indices, column_indices, trues(length(row_indices)),
                             n_state_dofs, n_state_dofs, |)
    if nnz(spatial_pattern) != expected_nonzeros
        throw(ArgumentError("Spatial Jacobian pattern contains duplicate or " *
                            "inconsistent entries."))
    end

    return spatial_pattern
end

# Retain the spatial pattern for the standard physical system
function physical_jacobian_sparsity_pattern(spatial_pattern,
                                            ::TemporalOperatorStandard)
    return spatial_pattern
end

# Add the nodal capacity derivative to the physical Jacobian pattern
function physical_jacobian_sparsity_pattern(spatial_pattern,
                                            ::TemporalOperatorCapacity)
    n_state_dofs = size(spatial_pattern, 1)
    identity_pattern = sparse(1:n_state_dofs, 1:n_state_dofs, trues(n_state_dofs),
                              n_state_dofs, n_state_dofs)
    return spatial_pattern .| identity_pattern
end

# Embed the spatial pattern in the constitutive evolved-state block layout
function physical_jacobian_sparsity_pattern(spatial_pattern,
                                            ::TemporalOperatorConstitutive)
    n_state_dofs = size(spatial_pattern, 1)
    zero_pattern = spzeros(Bool, n_state_dofs, n_state_dofs)
    identity_pattern = sparse(1:n_state_dofs, 1:n_state_dofs, trues(n_state_dofs),
                              n_state_dofs, n_state_dofs)
    return [zero_pattern spatial_pattern; identity_pattern identity_pattern]
end

# Convert the physical Jacobian pattern to the numerical prototype
function jacobian_prototype(u_ode, spatial_pattern, physical_pattern, solver, cache,
                            ::NoPassiveVariables)
    jac_prototype = SparseMatrixCSC{eltype(u_ode), Int}(physical_pattern)
    fill!(nonzeros(jac_prototype), zero(eltype(jac_prototype)))
    return jac_prototype
end

# Append boundary-flux rows that depend on the adjacent state elements
function jacobian_prototype(u_ode, spatial_pattern, physical_pattern, solver, cache,
                            ::PassiveVariablesBoundaryFlux1D)
    n_nodes = Trixi.nnodes(solver)
    n_state_dofs = size(spatial_pattern, 1)
    n_physical_dofs = size(physical_pattern, 1)
    state_offset = n_physical_dofs - n_state_dofs
    state_indices = LinearIndices((n_nodes, n_state_dofs ÷ n_nodes))

    n_boundaries_per_direction = cache.boundaries.n_boundaries_per_direction
    lasts = accumulate(+, n_boundaries_per_direction)
    firsts = lasts - n_boundaries_per_direction .+ 1
    element_neg = cache.boundaries.neighbor_ids[firsts[1]]
    element_pos = cache.boundaries.neighbor_ids[firsts[2]]

    boundary_pattern = spzeros(Bool, 2, n_physical_dofs)
    for state_node in Trixi.eachnode(solver)
        boundary_pattern[1, state_offset + state_indices[state_node, element_neg]] = true
        boundary_pattern[2, state_offset + state_indices[state_node, element_pos]] = true
    end

    passive_columns = spzeros(Bool, n_physical_dofs + 2, 2)
    augmented_pattern = [physical_pattern; boundary_pattern]
    augmented_pattern = [augmented_pattern passive_columns]
    jac_prototype = SparseMatrixCSC{eltype(u_ode), Int}(augmented_pattern)
    fill!(nonzeros(jac_prototype), zero(eltype(jac_prototype)))
    return jac_prototype
end
