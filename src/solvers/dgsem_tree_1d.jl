@muladd begin
#! format: noindent

# New parabolic boundary container to allow boundary conditions to depend on `u_inner`,
# which in base Trixi.jl is overwritten by the parabolic fluxes to reuse the same container.
# This creates a new container, allowing both to be available for use in the boundary flux.
mutable struct ParabolicBoundaryContainer1D{uEltype <: Real}
    flux_values::Array{uEltype, 3}
    _flux_values::Vector{uEltype}

    function ParabolicBoundaryContainer1D{uEltype}(capacity::Integer,
                                                   n_variables::Integer) where {uEltype <:
                                                                                Real}
        flux_values = Array{uEltype, 3}(undef, 2, n_variables, capacity)
        _flux_values = Vector{uEltype}(undef, 2 * n_variables * capacity)
        return new(flux_values, _flux_values)
    end
end

# Creates the parabolic cache, including parabolic-specific boundary container
function Trixi.create_cache_parabolic(mesh::Trixi.TreeMesh{1},
                                      equations::Trixi.AbstractEquationsParabolic{1},
                                      dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                                      n_elements,
                                      ::Type{uEltype}) where {uEltype <: Real}
    parabolic_container = Trixi.init_parabolic_container_1d(Trixi.nvariables(equations),
                                                            Trixi.nnodes(dg),
                                                            n_elements,
                                                            uEltype)
    n_boundaries = Trixi.isperiodic(mesh, 1) ? 0 : 2
    parabolic_boundaries = ParabolicBoundaryContainer1D{uEltype}(n_boundaries,
                                                                 Trixi.nvariables(equations))

    return (; parabolic_container, parabolic_boundaries)
end

function Base.resize!(boundaries::ParabolicBoundaryContainer1D, equations, dg, cache)
    capacity = 2 * Trixi.nvariables(equations) * Trixi.nboundaries(cache.boundaries)
    resize!(boundaries._flux_values, capacity)
    boundaries.flux_values = unsafe_wrap(Array,
                                         pointer(boundaries._flux_values),
                                         (2, Trixi.nvariables(equations),
                                          Trixi.nboundaries(cache.boundaries)))
    return nothing
end

function Trixi.refine!(u_ode::AbstractVector,
                       adaptor,
                       mesh::Trixi.TreeMesh{1},
                       equations::Trixi.AbstractEquationsParabolic{1},
                       dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                       cache,
                       cache_parabolic,
                       elements_to_refine,
                       limiter!)
    Trixi.refine!(u_ode, adaptor, mesh, equations, dg, cache, elements_to_refine,
                  limiter!)

    (; parabolic_container, parabolic_boundaries) = cache_parabolic
    resize!(parabolic_container, equations, dg, cache)
    resize!(parabolic_boundaries, equations, dg, cache)

    return nothing
end

function Trixi.coarsen!(u_ode::AbstractVector,
                        adaptor,
                        mesh::Trixi.TreeMesh{1},
                        equations::Trixi.AbstractEquationsParabolic{1},
                        dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                        cache,
                        cache_parabolic,
                        elements_to_remove,
                        limiter!)
    Trixi.coarsen!(u_ode, adaptor, mesh, equations, dg, cache, elements_to_remove,
                   limiter!)

    (; parabolic_container, parabolic_boundaries) = cache_parabolic
    resize!(parabolic_container, equations, dg, cache)
    resize!(parabolic_boundaries, equations, dg, cache)

    return nothing
end

# Version of `rhs_parabolic!` that allows both the interior-side solution and flux to be 
# available at the boundary for use in the boundary condition
function Trixi.rhs_parabolic!(du, u, t, mesh::Trixi.TreeMesh{1},
                              equations_parabolic::Trixi.AbstractEquationsParabolic{1},
                              boundary_conditions, source_terms,
                              solver::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                              solver_parabolic, cache, cache_parabolic)
    (; parabolic_container, parabolic_boundaries) = cache_parabolic
    (; u_transformed, gradients, flux_parabolic) = parabolic_container

    Trixi.@trixi_timeit Trixi.timer() "transform variables" begin
        Trixi.transform_variables!(u_transformed, u, mesh, equations_parabolic, solver,
                                   cache)
    end

    Trixi.@trixi_timeit Trixi.timer() "calculate gradient" begin
        Trixi.calc_gradient!(gradients, u_transformed, t, mesh, equations_parabolic,
                             boundary_conditions, solver, solver_parabolic, cache)
    end

    Trixi.@trixi_timeit Trixi.timer() "calculate parabolic fluxes" begin
        Trixi.calc_parabolic_fluxes!(flux_parabolic, gradients, u_transformed, mesh,
                                     equations_parabolic, solver, cache)
    end

    Trixi.@trixi_timeit Trixi.timer() "reset ∂u/∂t" Trixi.set_zero!(du, solver, cache)

    Trixi.@trixi_timeit Trixi.timer() "volume integral" begin
        Trixi.calc_volume_integral!(du, flux_parabolic, mesh, equations_parabolic,
                                    solver, cache)
    end

    Trixi.@trixi_timeit Trixi.timer() "prolong2interfaces" begin
        Trixi.prolong2interfaces!(cache, flux_parabolic, mesh, equations_parabolic,
                                  solver)
    end

    Trixi.@trixi_timeit Trixi.timer() "interface flux" begin
        Trixi.calc_interface_flux!(cache.elements.surface_flux_values,
                                   mesh, equations_parabolic, solver,
                                   solver_parabolic, cache)
    end

    # prolong2boundaries and boundary flux differ from Trixi.jl to allow both the 
    # interior-side solution and flux to be available at the boundary for use in the 
    # boundary condition
    Trixi.@trixi_timeit Trixi.timer() "prolong2boundaries" begin
        prolong2boundaries_divergence!(parabolic_boundaries, cache, flux_parabolic,
                                       mesh,
                                       equations_parabolic, solver)
    end

    Trixi.@trixi_timeit Trixi.timer() "boundary flux" begin
        calc_boundary_flux_divergence!(cache, t, boundary_conditions, mesh,
                                       equations_parabolic,
                                       solver.surface_integral, solver, cache_parabolic)
    end

    Trixi.@trixi_timeit Trixi.timer() "surface integral" begin
        Trixi.calc_surface_integral!(nothing, du, u, mesh, equations_parabolic,
                                     solver.surface_integral, solver, cache)
    end

    Trixi.@trixi_timeit Trixi.timer() "Jacobian" begin
        Trixi.apply_jacobian_parabolic!(du, mesh, equations_parabolic, solver, cache)
    end

    Trixi.@trixi_timeit Trixi.timer() "source terms parabolic" begin
        Trixi.calc_sources_parabolic!(du, u, gradients, t, source_terms,
                                      equations_parabolic, solver, cache)
    end

    return nothing
end

# If periodic, skip this boundary flux calculation, as there are no boundaries
function calc_boundary_flux_divergence!(cache, t,
                                        boundary_conditions::Trixi.BoundaryConditionPeriodic,
                                        mesh::Trixi.TreeMesh{1},
                                        equations_parabolic::Trixi.AbstractEquationsParabolic{1},
                                        surface_integral,
                                        dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                                        cache_parabolic)
    return nothing
end

function calc_boundary_flux_divergence!(cache, t,
                                        boundary_conditions::NamedTuple,
                                        mesh::Trixi.TreeMesh{1},
                                        equations_parabolic::Trixi.AbstractEquationsParabolic{1},
                                        surface_integral,
                                        dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                                        cache_parabolic)
    surface_flux_values = cache.elements.surface_flux_values
    n_boundaries_per_direction = cache.boundaries.n_boundaries_per_direction

    lasts = accumulate(+, n_boundaries_per_direction)
    firsts = lasts - n_boundaries_per_direction .+ 1

    calc_boundary_flux_by_direction_divergence!(surface_flux_values, t,
                                                boundary_conditions[1],
                                                equations_parabolic,
                                                surface_integral, dg,
                                                cache, 1, firsts[1], lasts[1],
                                                cache_parabolic)
    calc_boundary_flux_by_direction_divergence!(surface_flux_values, t,
                                                boundary_conditions[2],
                                                equations_parabolic,
                                                surface_integral, dg,
                                                cache, 2, firsts[2], lasts[2],
                                                cache_parabolic)

    return nothing
end

function calc_boundary_flux_by_direction_divergence!(surface_flux_values::AbstractArray{<:Any,
                                                                                        3},
                                                     t,
                                                     boundary_condition,
                                                     equations_parabolic::Trixi.AbstractEquationsParabolic{1},
                                                     surface_integral,
                                                     dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis},
                                                     cache, direction,
                                                     first_boundary, last_boundary,
                                                     cache_parabolic)
    parabolic_boundaries = cache_parabolic.parabolic_boundaries
    (; u, neighbor_ids, neighbor_sides, node_coordinates, orientations) = cache.boundaries
    parabolic_flux_values = parabolic_boundaries.flux_values

    Trixi.@threaded for boundary in first_boundary:last_boundary
        neighbor = neighbor_ids[boundary]

        flux_ll, flux_rr = Trixi.get_surface_node_vars(parabolic_flux_values,
                                                       equations_parabolic, dg,
                                                       boundary)
        u_ll, u_rr = Trixi.get_surface_node_vars(u, equations_parabolic, dg, boundary)
        if neighbor_sides[boundary] == 1
            flux_inner = flux_ll
            u_inner = u_ll
        else
            flux_inner = flux_rr
            u_inner = u_rr
        end

        x = Trixi.get_node_coords(node_coordinates, equations_parabolic, dg, boundary)
        flux = boundary_condition(flux_inner, u_inner, orientations[boundary],
                                  direction,
                                  x, t, Trixi.Divergence(), equations_parabolic)

        for v in Trixi.eachvariable(equations_parabolic)
            surface_flux_values[v, direction, neighbor] = flux[v]
        end
    end

    return nothing
end

function prolong2boundaries_divergence!(parabolic_boundaries, cache, flux_parabolic,
                                        mesh::Trixi.TreeMesh{1},
                                        equations::Trixi.AbstractEquationsParabolic{1},
                                        dg::Trixi.DGSEM{<:Trixi.LobattoLegendreBasis})
    neighbor_sides = cache.boundaries.neighbor_sides
    boundary_flux_values = parabolic_boundaries.flux_values

    Trixi.@threaded for boundary in Trixi.eachboundary(dg, cache)
        element = cache.boundaries.neighbor_ids[boundary]

        if neighbor_sides[boundary] == 1
            for v in Trixi.eachvariable(equations)
                boundary_flux_values[1, v, boundary] = flux_parabolic[v,
                                                                      Trixi.nnodes(dg),
                                                                      element]
            end
        else
            for v in Trixi.eachvariable(equations)
                boundary_flux_values[2, v, boundary] = flux_parabolic[v, 1, element]
            end
        end
    end

    return nothing
end
end # @muladd
