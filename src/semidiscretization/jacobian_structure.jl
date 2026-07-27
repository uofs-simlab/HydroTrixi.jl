function jacobian_structure(u_ode, semi::SemidiscretizationImplicit)
    semi_base = semi.semi_base
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi_base)

    return jacobian_structure(u_ode, mesh, equations, solver,
                              semi_base.solver_parabolic, cache,
                              semi.operator_temporal, semi.passive_variables)
end

function jacobian_structure(u_ode, mesh, equations, solver, solver_parabolic,
                            cache, operator_temporal, passive_variables)
    return nothing
end

function jacobian_structure(u_ode, mesh::Trixi.TreeMesh{1},
                            equations::RichardsEquation1D,
                            solver::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                            solver_parabolic::Trixi.ParabolicFormulationLocalDG{Nothing},
                            cache,
                            operator_temporal::TemporalOperatorConstitutive,
                            passive_variables::NoPassiveVariables)
    Trixi.mpi_isparallel() && return nothing
    Trixi.isperiodic(mesh) && return nothing

    n_nodes = Trixi.nnodes(solver)
    n_elements = Trixi.nelements(solver, cache)
    n_interfaces = length(Trixi.eachinterface(solver, cache))
    n_physical_dofs = n_nodes * n_elements
    n_total_dofs = 2 * n_physical_dofs
    length(u_ode) == n_total_dofs ||
        throw(DimensionMismatch("Expected $n_total_dofs mixed-form degrees of freedom, " *
                                "received $(length(u_ode))."))

    expected_nonzeros = n_elements * n_nodes^2 + 2 * n_nodes * n_interfaces +
                        2 * n_physical_dofs
    row_indices = Int[]
    column_indices = Int[]
    sizehint!(row_indices, expected_nonzeros)
    sizehint!(column_indices, expected_nonzeros)
    physical_indices = LinearIndices((n_nodes, n_elements))

    # Add the mixed identity and constitutive block diagonals
    for physical_dof in 1:n_physical_dofs
        push!(row_indices, n_physical_dofs + physical_dof)
        push!(column_indices, physical_dof)
        push!(row_indices, n_physical_dofs + physical_dof)
        push!(column_indices, n_physical_dofs + physical_dof)
    end

    # Add dense element-local spatial dependencies
    for element in Trixi.eachelement(solver, cache)
        for state_node in Trixi.eachnode(solver)
            state_dof = physical_indices[state_node, element]
            for residual_node in Trixi.eachnode(solver)
                residual_dof = physical_indices[residual_node, element]
                push!(row_indices, residual_dof)
                push!(column_indices, n_physical_dofs + state_dof)
            end
        end
    end

    # Add the asymmetric LDG interface dependencies
    for interface in Trixi.eachinterface(solver, cache)
        left_element = cache.interfaces.neighbor_ids[1, interface]
        right_element = cache.interfaces.neighbor_ids[2, interface]
        left_face_dof = physical_indices[n_nodes, left_element]

        for right_node in Trixi.eachnode(solver)
            right_dof = physical_indices[right_node, right_element]
            push!(row_indices, left_face_dof)
            push!(column_indices, n_physical_dofs + right_dof)
            push!(row_indices, right_dof)
            push!(column_indices, n_physical_dofs + left_face_dof)
        end
    end

    jac_prototype = sparse(row_indices, column_indices,
                           zeros(eltype(u_ode), length(row_indices)),
                           n_total_dofs, n_total_dofs)
    nnz(jac_prototype) == expected_nonzeros ||
        throw(ArgumentError("Sparse Jacobian structure contains duplicate entries."))

    return jac_prototype
end
