@inline _uses_penalized_dirichlet_boundary_condition(::Any) = false
@inline _uses_penalized_dirichlet_boundary_condition(::BoundaryConditionDirichletPenalty) = true
@inline _uses_penalized_dirichlet_boundary_condition(boundary_conditions::NamedTuple) =
    any(_uses_penalized_dirichlet_boundary_condition, values(boundary_conditions))

function rhs_parabolic!(du, u, t, mesh, equations, boundary_conditions, source_terms,
                        solver, parabolic_scheme, base_cache, viscous_cache)
    if _uses_penalized_dirichlet_boundary_condition(boundary_conditions)
        if mesh isa Trixi.TreeMesh{1}
            rhs_parabolic_penalized_ldg!(du, u, t, mesh, equations, boundary_conditions,
                                         source_terms, solver, parabolic_scheme, base_cache,
                                         viscous_cache)
        else
            throw(ArgumentError("BoundaryConditionDirichletPenalty currently supports " *
                                "TreeMesh{1} semidiscretizations."))
        end
    else
        Trixi.rhs_parabolic!(du, u, t, mesh, equations, boundary_conditions, source_terms,
                             solver, parabolic_scheme, base_cache, viscous_cache)
    end

    return nothing
end

function rhs_parabolic_penalized_ldg!(du, u, t,
                                      mesh::Trixi.TreeMesh{1},
                                      equations::Trixi.AbstractEquationsParabolic,
                                      boundary_conditions::NamedTuple,
                                      source_terms,
                                      solver::Trixi.DG,
                                      parabolic_scheme,
                                      base_cache, viscous_cache)
    viscous_container = viscous_cache.viscous_container
    u_transformed = viscous_container.u_transformed
    gradients = viscous_container.gradients
    flux_viscous = viscous_container.flux_viscous

    Trixi.transform_variables!(u_transformed, u, mesh, equations, solver, base_cache)
    Trixi.calc_gradient!(gradients, u_transformed, t, mesh, equations, boundary_conditions,
                         solver, parabolic_scheme, base_cache)

    # Preserve boundary traces of the solution before `cache.boundaries.u` is overwritten
    # with viscous fluxes in the divergence stage.
    boundary_solution_values = copy(base_cache.boundaries.u)

    Trixi.calc_viscous_fluxes!(flux_viscous, gradients, u_transformed, mesh,
                               equations, solver, base_cache)
    du .= zero(eltype(du))
    Trixi.calc_volume_integral!(du, flux_viscous, mesh, equations, solver, base_cache)
    Trixi.prolong2interfaces!(base_cache, flux_viscous, mesh, equations, solver)
    Trixi.calc_interface_flux!(base_cache.elements.surface_flux_values,
                               mesh, equations, solver, parabolic_scheme, base_cache)
    Trixi.prolong2boundaries!(base_cache, flux_viscous, mesh, equations, solver)

    calc_boundary_flux_divergence_with_solution!(base_cache, boundary_solution_values, t,
                                                 boundary_conditions, mesh, equations,
                                                 solver.surface_integral, solver)

    Trixi.calc_surface_integral!(du, u, mesh, equations,
                                 solver.surface_integral, solver, base_cache)
    Trixi.apply_jacobian_parabolic!(du, mesh, equations, solver, base_cache)
    Trixi.calc_sources_parabolic!(du, u, gradients, t, source_terms,
                                  equations, solver, base_cache)

    return nothing
end

function calc_boundary_flux_divergence_with_solution!(cache, boundary_solution_values, t,
                                                      boundary_conditions::NamedTuple,
                                                      mesh::Trixi.TreeMesh{1},
                                                      equations::Trixi.AbstractEquationsParabolic,
                                                      surface_integral, solver::Trixi.DG)
    surface_flux_values = cache.elements.surface_flux_values
    n_boundaries_per_direction = cache.boundaries.n_boundaries_per_direction

    lasts = accumulate(+, n_boundaries_per_direction)
    firsts = lasts - n_boundaries_per_direction .+ 1

    calc_boundary_flux_by_direction_divergence_with_solution!(surface_flux_values,
                                                              boundary_solution_values, t,
                                                              boundary_conditions[1], equations,
                                                              surface_integral, solver, cache,
                                                              1, firsts[1], lasts[1])
    calc_boundary_flux_by_direction_divergence_with_solution!(surface_flux_values,
                                                              boundary_solution_values, t,
                                                              boundary_conditions[2], equations,
                                                              surface_integral, solver, cache,
                                                              2, firsts[2], lasts[2])

    return nothing
end

function calc_boundary_flux_by_direction_divergence_with_solution!(surface_flux_values,
                                                                   boundary_solution_values, t,
                                                                   boundary_condition,
                                                                   equations::Trixi.AbstractEquationsParabolic,
                                                                   surface_integral,
                                                                   solver::Trixi.DG, cache,
                                                                   direction, first_boundary,
                                                                   last_boundary)
    flux_boundary_values = cache.boundaries.u
    neighbor_ids = cache.boundaries.neighbor_ids
    neighbor_sides = cache.boundaries.neighbor_sides
    node_coordinates = cache.boundaries.node_coordinates
    orientations = cache.boundaries.orientations

    for boundary in first_boundary:last_boundary
        neighbor = neighbor_ids[boundary]

        flux_ll, flux_rr = Trixi.get_surface_node_vars(flux_boundary_values, equations, solver,
                                                       boundary)
        u_ll, u_rr = Trixi.get_surface_node_vars(boundary_solution_values, equations, solver,
                                                 boundary)

        if neighbor_sides[boundary] == 1
            flux_inner = flux_ll
            u_inner = u_ll
        else
            flux_inner = flux_rr
            u_inner = u_rr
        end

        x = Trixi.get_node_coords(node_coordinates, equations, solver, boundary)
        flux = boundary_condition(flux_inner, u_inner, orientations[boundary], direction,
                                  x, t, Trixi.Divergence(), equations)

        for v in Trixi.eachvariable(equations)
            surface_flux_values[v, direction, neighbor] = flux[v]
        end
    end

    return nothing
end
