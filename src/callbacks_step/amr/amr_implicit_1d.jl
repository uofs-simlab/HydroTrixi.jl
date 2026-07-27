function fill_mass_matrix_diagonal!(diagonal_entries,
                                    ::TemporalOperatorConstitutive)
    fill!(diagonal_entries, zero(eltype(diagonal_entries)))

    half = length(diagonal_entries) ÷ 2
    @inbounds diagonal_entries[1:half] .= one(eltype(diagonal_entries))

    return nothing
end

function fill_mass_matrix_diagonal!(diagonal_entries, ::AbstractTemporalOperator)
    fill!(diagonal_entries, one(eltype(diagonal_entries)))
    return nothing
end

function update_mass_matrix!(mass_matrix::Diagonal, u_ode,
                             operator_temporal::TemporalOperatorConstitutive,
                             semi_base::Trixi.AbstractSemidiscretization)
    diagonal_entries = mass_matrix.diag
    resize!(diagonal_entries, length(u_ode))
    fill_mass_matrix_diagonal!(diagonal_entries, operator_temporal)

    return nothing
end

function update_mass_matrix!(mass_matrix_implicit::Diagonal, u_ode,
                             semi::SemidiscretizationImplicit)
    n_passive = passive_variable_count(semi)
    n_physical = length(u_ode) - n_passive
    diagonal_entries = mass_matrix_implicit.diag
    resize!(diagonal_entries, length(u_ode))

    # The physical block keeps the temporal-operator mass matrix layout
    diagonal_entries_physical = @view diagonal_entries[1:n_physical]
    fill_mass_matrix_diagonal!(diagonal_entries_physical, semi.operator_temporal)

    if n_passive > 0
        # Passive scalar variables are differential variables
        diagonal_entries_passive = @view diagonal_entries[(n_physical + 1):end]
        fill!(diagonal_entries_passive, one(eltype(diagonal_entries)))
    end

    return nothing
end

function update_mass_matrix!(mass_matrix, u_ode, operator_temporal,
                             semi_base::Trixi.AbstractSemidiscretization)
    return nothing
end

function update_mass_matrix!(mass_matrix, u_ode, semi::SemidiscretizationImplicit)
    return nothing
end

function update_mass_matrix!(ode_function::SciMLBase.ODEFunction, u_ode,
                             operator_temporal::AbstractTemporalOperator,
                             semi_base::Trixi.AbstractSemidiscretization)
    update_mass_matrix!(ode_function.mass_matrix, u_ode, operator_temporal, semi_base)

    if ode_function.f isa SciMLBase.ODEFunction
        update_mass_matrix!(ode_function.f, u_ode, operator_temporal, semi_base)
    end

    return nothing
end

function update_mass_matrix!(ode_function::SciMLBase.ODEFunction, u_ode,
                             semi::SemidiscretizationImplicit)
    update_mass_matrix!(ode_function.mass_matrix, u_ode, semi)

    if ode_function.f isa SciMLBase.ODEFunction
        update_mass_matrix!(ode_function.f, u_ode, semi)
    end

    return nothing
end

# Rebuild factorization caches after AMR changes the state size
function refresh_linear_solver_cache!(cache)
    if !(hasproperty(cache, :linsolve) && hasproperty(cache, :W) &&
         hasproperty(cache, :linsolve_tmp) && hasproperty(cache, :tmp))
        return nothing
    end

    linsolve = cache.linsolve
    linear_problem = SciMLBase.LinearProblem(cache.W, vec(cache.linsolve_tmp),
                                             linsolve.p; u0 = vec(cache.tmp))
    alias = SciMLBase.LinearAliasSpecifier(alias_A = true, alias_b = true)
    cache.linsolve = SciMLBase.init(linear_problem, linsolve.alg;
                                    alias = alias,
                                    Pl = linsolve.Pl,
                                    Pr = linsolve.Pr,
                                    abstol = linsolve.abstol,
                                    reltol = linsolve.reltol,
                                    maxiters = linsolve.maxiters,
                                    verbose = linsolve.verbose,
                                    assumptions = linsolve.assumptions,
                                    sensealg = linsolve.sensealg)

    return nothing
end

struct AMRCallbackImplicit{Controller, Adaptor, Cache}
    controller::Controller
    interval::Int
    adapt_initial_condition::Bool
    adapt_initial_condition_only_refine::Bool
    dynamic_load_balancing::Bool
    adaptor::Adaptor
    amr_cache::Cache
end

function Trixi.AMRCallback(semi::SemidiscretizationImplicit, controller, adaptor;
                           interval,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true,
                           dynamic_load_balancing = true)
    semi.operator_temporal isa TemporalOperatorConstitutive &&
        isnothing(semi.operator_temporal.evolved_to_state) &&
        throw(ArgumentError("Constitutive AMR requires `evolved_to_state`."))

    if !(interval isa Integer && interval >= 0)
        throw(ArgumentError("`interval` must be a non-negative integer."))
    end

    if interval > 0
        condition = (u, t, integrator) -> ((integrator.stats.naccept % interval == 0) &&
                                           !(integrator.stats.naccept == 0 &&
                                             integrator.iter > 0) &&
                                           !Trixi.isfinished(integrator))
    else
        condition = (u, t, integrator) -> false
    end

    amr_cache = (; to_refine = Int[], to_coarsen = Int[])
    amr_callback = AMRCallbackImplicit(controller,
                                       Int(interval),
                                       adapt_initial_condition,
                                       adapt_initial_condition_only_refine,
                                       dynamic_load_balancing,
                                       adaptor,
                                       amr_cache)

    return SciMLBase.DiscreteCallback(condition, amr_callback;
                                      save_positions = (false, false),
                                      initialize = initialize_amr!)
end

function Trixi.AMRCallback(semi::SemidiscretizationImplicit, controller; kwargs...)
    adaptor = Trixi.AdaptorAMR(semi.semi_base)
    return Trixi.AMRCallback(semi, controller, adaptor; kwargs...)
end

function initialize_amr!(cb::SciMLBase.DiscreteCallback{Condition, Affect!}, u, t,
                         integrator) where {Condition,
                                            Affect! <: AMRCallbackImplicit}
    amr_callback = cb.affect!
    semi = integrator.p

    if amr_callback.adapt_initial_condition
        only_refine = amr_callback.adapt_initial_condition_only_refine
        has_changed = amr_callback(integrator;
                                   only_refine = only_refine)
        iterations = 1
        while has_changed
            Trixi.compute_coefficients!(integrator.u, t, semi)
            SciMLBase.u_modified!(integrator, true)
            has_changed = amr_callback(integrator;
                                       only_refine = only_refine)
            iterations += 1
            allowed_max_iterations = max(10, Trixi.max_level(amr_callback.controller))
            if iterations > allowed_max_iterations
                @warn "AMR for initial condition did not settle within " *
                      "$(allowed_max_iterations) iterations."
                break
            end
        end
    end

    return nothing
end

function (amr_callback::AMRCallbackImplicit)(integrator; kwargs...)
    isnothing(integrator.f.jac_prototype) ||
        throw(ArgumentError("Sparse Jacobian structure cannot be used with adaptive " *
                            "mesh refinement until the structure can be rebuilt."))

    u_ode = integrator.u
    semi = integrator.p

    has_changed = amr_callback(u_ode, semi, integrator.t, integrator.iter; kwargs...)
    if has_changed
        update_mass_matrix!(integrator.f, u_ode, semi)
        resize!(integrator, length(u_ode))
        refresh_linear_solver_cache!(integrator.cache)
        SciMLBase.u_modified!(integrator, true)
    end

    return has_changed
end

function (amr_callback::AMRCallbackImplicit)(u_ode::AbstractVector,
                                             semi::SemidiscretizationImplicit,
                                             t, iter;
                                             only_refine = false,
                                             only_coarsen = false)
    semi_base = semi.semi_base
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi_base)

    if Trixi.mpi_isparallel()
        error("MPI AMR has not been verified for `SemidiscretizationImplicit`.")
    end

    u_state = Trixi.wrap_array(state_variable_view(u_ode, semi),
                               mesh, equations, dg, cache)
    lambda = amr_callback.controller(u_state, mesh, equations, dg, cache;
                                     t = t, iter = iter)
    leaf_cell_ids = Trixi.leaf_cells(mesh.tree)

    @boundscheck begin
        @assert axes(lambda)==axes(leaf_cell_ids) ("Indicator and leaf cell arrays "*
                                                   "have different axes")
    end

    (; to_refine, to_coarsen) = amr_callback.amr_cache
    empty!(to_refine)
    empty!(to_coarsen)

    for element in eachindex(lambda)
        controller_value = lambda[element]
        if controller_value > 0
            push!(to_refine, leaf_cell_ids[element])
        elseif controller_value < 0
            push!(to_coarsen, leaf_cell_ids[element])
        end
    end

    transferred_ode = transferred_variables_for_amr(u_ode, semi,
                                                    semi.operator_temporal)
    passive_ode = collect(passive_variable_view(u_ode, semi))
    refined_original_cells = refine_transferred_variables!(transferred_ode,
                                                           amr_callback, semi_base,
                                                           to_refine, only_coarsen)
    coarsened_original_cells = coarsen_transferred_variables!(transferred_ode,
                                                              amr_callback, semi_base,
                                                              to_coarsen,
                                                              refined_original_cells,
                                                              only_refine)

    has_changed = !isempty(refined_original_cells) || !isempty(coarsened_original_cells)
    if has_changed
        resize_after_amr!(u_ode, transferred_ode, passive_ode, semi,
                          semi.operator_temporal)
        mesh.unsaved_changes = true
    end

    return has_changed
end

function transferred_variables_for_amr(u_ode, semi::SemidiscretizationImplicit,
                                       ::AbstractTemporalOperator)
    return collect(evolved_variable_view(u_ode, semi))
end

function transferred_variables_for_amr(u_ode, semi::SemidiscretizationImplicit,
                                       operator_temporal::TemporalOperatorCapacity)
    transferred_ode = collect(state_variable_view(u_ode, semi))
    equations = semi.semi_base.equations
    transfer_variables = operator_temporal.transfer_variables

    # Convert the state to the variable used for AMR grid transfer
    @inbounds for i in eachindex(transferred_ode)
        transferred_ode[i] = transfer_variables(transferred_ode[i], equations)
    end

    return transferred_ode
end

function resize_after_amr!(u_ode, transferred_ode, passive_ode,
                           semi::SemidiscretizationImplicit,
                           operator_temporal::TemporalOperatorCapacity)
    # Passive scalar variables are global diagnostics and are not adapted
    resize!(u_ode, length(transferred_ode) + length(passive_ode))
    state_variable = physical_variable_view(u_ode, semi)
    equations = semi.semi_base.equations
    transfer_to_state = operator_temporal.transfer_to_state

    # Reconstruct the pressure head state from the transferred variable
    @inbounds for i in eachindex(state_variable, transferred_ode)
        state_variable[i] = transfer_to_state(transferred_ode[i], equations)
    end

    passive_variable_view(u_ode, semi) .= passive_ode
    return nothing
end

function resize_after_amr!(u_ode, transferred_ode, passive_ode,
                           semi::SemidiscretizationImplicit,
                           ::TemporalOperatorConstitutive)
    # Passive scalar variables are global diagnostics and are not adapted
    resize!(u_ode, 2 * length(transferred_ode) + length(passive_ode))
    evolved_variable_view(u_ode, semi) .= transferred_ode
    reconstruct_state_from_evolved!(u_ode, semi)
    passive_variable_view(u_ode, semi) .= passive_ode
    return nothing
end

function refine_transferred_variables!(transferred_ode,
                                       amr_callback::AMRCallbackImplicit,
                                       semi_base::Trixi.AbstractSemidiscretization,
                                       to_refine, only_coarsen)
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi_base)

    if only_coarsen || isempty(to_refine)
        return Int[]
    end

    refined_original_cells = Trixi.refine!(mesh.tree, to_refine)
    elements_to_refine = findall(in(refined_original_cells), cache.elements.cell_ids)
    Trixi.refine!(transferred_ode, amr_callback.adaptor, mesh, equations, dg, cache,
                  semi_base.cache_parabolic, elements_to_refine)

    return refined_original_cells
end

function coarsen_transferred_variables!(transferred_ode,
                                        amr_callback::AMRCallbackImplicit,
                                        semi_base::Trixi.AbstractSemidiscretization,
                                        to_coarsen, refined_original_cells, only_refine)
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi_base)

    if only_refine || isempty(to_coarsen)
        return Int[]
    end

    if !isempty(to_coarsen)
        to_coarsen = Trixi.original2refined(to_coarsen, refined_original_cells, mesh)
    end

    parents_to_coarsen = zeros(Int, length(mesh.tree))
    for cell_id in to_coarsen
        if !Trixi.has_parent(mesh.tree, cell_id) || !Trixi.is_leaf(mesh.tree, cell_id)
            continue
        end

        parent_id = mesh.tree.parent_ids[cell_id]
        parents_to_coarsen[parent_id] += 1
    end

    n_children = Trixi.n_children_per_cell(mesh.tree)
    to_coarsen = collect(eachindex(parents_to_coarsen))[parents_to_coarsen .== n_children]
    coarsened_original_cells = Trixi.coarsen!(mesh.tree, to_coarsen)

    removed_child_cells = zeros(Int, n_children * length(coarsened_original_cells))
    for (index, coarse_cell_id) in enumerate(coarsened_original_cells)
        for child in 1:n_children
            removed_child_cells[n_children * (index - 1) + child] = coarse_cell_id + child
        end
    end

    elements_to_remove = findall(in(removed_child_cells), cache.elements.cell_ids)
    Trixi.coarsen!(transferred_ode, amr_callback.adaptor, mesh, equations, dg, cache,
                   semi_base.cache_parabolic, elements_to_remove)

    return coarsened_original_cells
end

function reconstruct_state_from_evolved!(u_ode, semi::SemidiscretizationImplicit)
    evolved_variable = evolved_variable_view(u_ode, semi)
    state_variable = state_variable_view(u_ode, semi)
    equations = semi.semi_base.equations
    evolved_to_state = semi.operator_temporal.evolved_to_state
    isnothing(evolved_to_state) &&
        throw(ArgumentError("Constitutive AMR requires `evolved_to_state`."))

    @inbounds for i in eachindex(state_variable, evolved_variable)
        state_variable[i] = evolved_to_state(evolved_variable[i], equations)
    end

    return nothing
end
