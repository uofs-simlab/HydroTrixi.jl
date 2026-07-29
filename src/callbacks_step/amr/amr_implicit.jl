# Mark evolved variables as differential and state variables as algebraic
function fill_mass_matrix_diagonal!(diagonal_entries, ::TemporalOperatorConstitutive)
    fill!(diagonal_entries, zero(eltype(diagonal_entries)))
    half = length(diagonal_entries) ÷ 2
    @inbounds diagonal_entries[1:half] .= one(eltype(diagonal_entries))
    return nothing
end

# Mark all variables as differential for standard temporal operators
function fill_mass_matrix_diagonal!(diagonal_entries, ::AbstractTemporalOperator)
    fill!(diagonal_entries, one(eltype(diagonal_entries)))
    return nothing
end

# Resize the physical and passive blocks of an implicit mass matrix after AMR
function update_mass_matrix!(mass_matrix_implicit::Diagonal, u_ode,
                             semi::SemidiscretizationImplicit)
    n_passive = passive_variable_count(semi)
    n_physical = length(u_ode) - n_passive
    resize!(mass_matrix_implicit.diag, length(u_ode))

    # The physical block keeps the temporal-operator mass matrix layout
    diagonal_entries_physical = @view mass_matrix_implicit.diag[1:n_physical]
    fill_mass_matrix_diagonal!(diagonal_entries_physical, semi.operator_temporal)

    if n_passive > 0
        # Passive scalar variables are differential variables
        diagonal_entries_passive = @view mass_matrix_implicit.diag[(n_physical + 1):end]
        fill!(diagonal_entries_passive, one(eltype(mass_matrix_implicit.diag)))
    end

    return nothing
end

# Leave implicit mass matrices without resizable diagonal storage unchanged
function update_mass_matrix!(mass_matrix, u_ode, semi::SemidiscretizationImplicit)
    return nothing
end

# Propagate implicit mass-matrix updates through nested ODE function wrappers
function update_mass_matrix!(ode_function::SciMLBase.ODEFunction, u_ode,
                             semi::SemidiscretizationImplicit)
    update_mass_matrix!(ode_function.mass_matrix, u_ode, semi)

    if ode_function.f isa SciMLBase.ODEFunction
        update_mass_matrix!(ode_function.f, u_ode, semi)
    end

    return nothing
end

# Leave solver caches without a reusable linear solve unchanged
refresh_linear_solver_cache!(cache) = nothing

# Rebuild the Rosenbrock factorization cache after AMR changes the state size
function refresh_linear_solver_cache!(cache::OrdinaryDiffEqRosenbrock.RosenbrockCache)
    linsolve = cache.linsolve
    linear_problem = SciMLBase.LinearProblem(cache.W, vec(cache.linsolve_tmp),
                                             linsolve.p; u0 = vec(cache.tmp))
    cache.linsolve = SciMLBase.init(linear_problem, linsolve.alg;
                                    alias = SciMLBase.LinearAliasSpecifier(alias_A = true,
                                                                           alias_b = true),
                                    Pl = linsolve.Pl, Pr = linsolve.Pr,
                                    abstol = linsolve.abstol, reltol = linsolve.reltol,
                                    maxiters = linsolve.maxiters,
                                    verbose = linsolve.verbose,
                                    assumptions = linsolve.assumptions,
                                    sensealg = linsolve.sensealg)

    return nothing
end

# Install the adapted sparsity pattern and force a new Rosenbrock Jacobian
function refresh_linear_solver_cache!(cache::OrdinaryDiffEqRosenbrock.RosenbrockCache,
                                      jac_prototype, dt)
    cache.J = copy(jac_prototype)
    fill!(nonzeros(cache.J), zero(eltype(cache.J)))
    cache.W = copy(cache.J)
    cache.jac_reuse = OrdinaryDiffEqRosenbrock.JacReuseState(zero(dt),
                                                             cache.alg.max_jac_age)

    return refresh_linear_solver_cache!(cache)
end

# Store the AMR policy, transfer adaptor, and reusable cell-index buffers
struct AMRCallbackImplicit{Controller, Adaptor, Cache}
    controller::Controller
    interval::Int
    adapt_initial_condition::Bool
    adapt_initial_condition_only_refine::Bool
    dynamic_load_balancing::Bool
    adaptor::Adaptor
    amr_cache::Cache
end

# Construct an implicit AMR callback with an explicitly configured transfer adaptor
function Trixi.AMRCallback(semi::SemidiscretizationImplicit, controller, adaptor; interval,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true,
                           dynamic_load_balancing = true)
    interval = Int(interval)
    if interval < 0
        throw(ArgumentError("`interval` must be non-negative."))
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
    amr_callback = AMRCallbackImplicit(controller, interval, adapt_initial_condition,
                                       adapt_initial_condition_only_refine,
                                       dynamic_load_balancing, adaptor, amr_cache)

    return SciMLBase.DiscreteCallback(condition, amr_callback;
                                      save_positions = (false, false),
                                      initialize = initialize_amr!)
end

# Construct an implicit AMR callback with the default adaptor for the base discretization
function Trixi.AMRCallback(semi::SemidiscretizationImplicit, controller; kwargs...)
    adaptor = Trixi.AdaptorAMR(semi.semi_base)
    return Trixi.AMRCallback(semi, controller, adaptor; kwargs...)
end

# Adapt the initial condition repeatedly until the controller accepts the mesh
function initialize_amr!(cb::SciMLBase.DiscreteCallback{Condition, Affect!}, u, t,
                         integrator) where {Condition, Affect! <: AMRCallbackImplicit}
    amr_callback = cb.affect!
    semi = integrator.p

    if amr_callback.adapt_initial_condition
        only_refine = amr_callback.adapt_initial_condition_only_refine
        has_changed = amr_callback(integrator; only_refine = only_refine)
        iterations = 1
        while has_changed
            Trixi.compute_coefficients!(integrator.u, t, semi)
            SciMLBase.u_modified!(integrator, true)
            has_changed = amr_callback(integrator; only_refine = only_refine)
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

# Apply AMR to an integrator and rebuild storage whose size depends on the mesh
function (amr_callback::AMRCallbackImplicit)(integrator; kwargs...)
    u_ode = integrator.u
    semi = integrator.p
    has_jac_prototype = !isnothing(integrator.f.jac_prototype)

    # If the mesh has changed, update the mass matrix and rebuild the ODE problem
    has_changed = amr_callback(u_ode, semi, integrator.t, integrator.iter; kwargs...)
    if has_changed
        update_mass_matrix!(integrator.f, u_ode, semi)
        if has_jac_prototype
            # Rebuild the augmented residual prototype for the adapted mesh
            residual_prototype = residual_jacobian_prototype(integrator.u, integrator.p)
            expected_size = (length(integrator.u), length(integrator.u))
            if size(residual_prototype) != expected_size ||
               !all(iszero, nonzeros(residual_prototype))
                throw(ArgumentError("Adapted residual Jacobian prototype is inconsistent " *
                                    "with the ODE state."))
            end
            ode_function = SciMLBase.remake(integrator.f;
                                            sparsity = copy(residual_prototype),
                                            jac_prototype = copy(residual_prototype))
            ode_problem = SciMLBase.remake(integrator.sol.prob;
                                           f = ode_function, u0 = copy(integrator.u))

            # Reuse OrdinaryDiffEq's preparation for the adapted sparsity pattern
            unprepared_algorithm = SciMLBase.remake(integrator.alg;
                                                    autodiff = integrator.alg.autodiff.dense_ad)
            algorithm = OrdinaryDiffEq.OrdinaryDiffEqCore.DiffEqBase.prepare_alg(unprepared_algorithm,
                                                                                 integrator.u,
                                                                                 integrator.p,
                                                                                 ode_problem)
            @assert typeof(algorithm)===typeof(integrator.alg) ("Sparse Jacobian AMR "*
                                                                "changed the ODE "*
                                                                "algorithm type")

            # Synchronize the solution metadata with the rebuilt problem
            ode_function = ode_problem.f
            interpolation = OrdinaryDiffEq.OrdinaryDiffEqCore.InterpolationData(integrator.sol.interp,
                                                                                ode_function)
            solution = integrator.sol
            solution = SciMLBase.@set solution.prob = ode_problem
            solution = SciMLBase.@set solution.alg = algorithm
            solution = SciMLBase.@set solution.interp = interpolation

            # Install the rebuilt objects in the integrator and stage functions
            integrator.f = ode_function
            integrator.alg = algorithm
            integrator.sol = solution
            integrator.cache.alg = algorithm
            integrator.cache.tf.f = ode_function
            integrator.cache.uf.f = ode_function

            resize!(integrator, length(integrator.u))
            refresh_linear_solver_cache!(integrator.cache, ode_function.jac_prototype,
                                         integrator.dt)
        else
            resize!(integrator, length(u_ode))
            refresh_linear_solver_cache!(integrator.cache)
        end
        SciMLBase.u_modified!(integrator, true)
    end

    return has_changed
end

# Select cells, transfer physical variables, and resize the ODE state after AMR
function (amr_callback::AMRCallbackImplicit)(u_ode::AbstractVector,
                                             semi::SemidiscretizationImplicit, t, iter;
                                             only_refine = false, only_coarsen = false)
    semi_base = semi.semi_base
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi_base)

    if Trixi.mpi_isparallel()
        error("MPI AMR has not been verified for `SemidiscretizationImplicit`.")
    end

    # Evaluate the controller on the physical state represented on the current mesh
    u_state = Trixi.wrap_array(state_variable_view(u_ode, semi), mesh, equations, dg, cache)
    lambda = amr_callback.controller(u_state, mesh, equations, dg, cache;
                                     t = t, iter = iter)
    leaf_cell_ids = Trixi.leaf_cells(mesh.tree)

    @boundscheck begin
        @assert axes(lambda)==axes(leaf_cell_ids) ("Indicator and leaf cell arrays "*
                                                   "have different axes")
    end

    # Translate signed controller values into tree cell IDs for adaptation
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

    # Transfer the operator-specific physical variable while preserving passive values
    transferred_ode = transferred_variables_for_amr(u_ode, semi, semi.operator_temporal)
    passive_ode = collect(passive_variable_view(u_ode, semi))

    # Refine requested cells and prolong the physical transfer variable
    if only_coarsen || isempty(to_refine)
        refined_original_cells = Int[]
    else
        refined_original_cells = Trixi.refine!(mesh.tree, to_refine)
        elements_to_refine = findall(in(refined_original_cells), cache.elements.cell_ids)
        Trixi.refine!(transferred_ode, amr_callback.adaptor, mesh, equations, dg,
                      cache, semi_base.cache_parabolic, elements_to_refine, nothing)
    end

    # Coarsen complete sibling groups and restrict the physical transfer variable
    if only_refine || isempty(to_coarsen)
        coarsened_original_cells = Int[]
    else
        mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi_base)
        to_coarsen = Trixi.original2refined(to_coarsen, refined_original_cells, mesh)

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
                removed_child_cells[n_children * (index - 1) + child] = coarse_cell_id +
                                                                        child
            end
        end

        elements_to_remove = findall(in(removed_child_cells), cache.elements.cell_ids)
        Trixi.coarsen!(transferred_ode, amr_callback.adaptor, mesh, equations, dg,
                       cache, semi_base.cache_parabolic, elements_to_remove, nothing)
    end

    has_changed = !isempty(refined_original_cells) || !isempty(coarsened_original_cells)
    if has_changed
        resize_after_amr!(u_ode, transferred_ode, passive_ode, semi, semi.operator_temporal)
        mesh.unsaved_changes = true
    end

    return has_changed
end

# Transfer evolved variables directly for standard and constitutive operators
function transferred_variables_for_amr(u_ode, semi::SemidiscretizationImplicit,
                                       ::Union{TemporalOperatorStandard,
                                               TemporalOperatorConstitutive})
    return collect(evolved_variable_view(u_ode, semi))
end

# Convert a capacity-form state to its configured AMR transfer variable
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

# Rebuild a capacity-form state from adapted transfer variables
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

# Rebuild the standard state while preserving passive diagnostic values
function resize_after_amr!(u_ode, transferred_ode, passive_ode,
                           semi::SemidiscretizationImplicit, ::TemporalOperatorStandard)
    resize!(u_ode, length(transferred_ode) + length(passive_ode))
    physical_variable_view(u_ode, semi) .= transferred_ode
    passive_variable_view(u_ode, semi) .= passive_ode
    return nothing
end

# Rebuild both constitutive blocks from the adapted evolved variables
function resize_after_amr!(u_ode, transferred_ode, passive_ode,
                           semi::SemidiscretizationImplicit, ::TemporalOperatorConstitutive)
    # Passive scalar variables are global diagnostics and are not adapted
    resize!(u_ode, 2 * length(transferred_ode) + length(passive_ode))
    evolved_variable = evolved_variable_view(u_ode, semi)
    evolved_variable .= transferred_ode
    state_variable = state_variable_view(u_ode, semi)
    equations = semi.semi_base.equations
    evolved_to_state = semi.operator_temporal.evolved_to_state

    @inbounds for i in eachindex(state_variable, evolved_variable)
        state_variable[i] = evolved_to_state(evolved_variable[i], equations)
    end

    passive_variable_view(u_ode, semi) .= passive_ode
    return nothing
end
