abstract type AbstractTemporalOperator end

@doc raw"""
    NoPassiveVariables()

Passive-variable configuration for implicit semidiscretizations without appended
diagnostic variables.
"""
struct NoPassiveVariables end

@doc raw"""
    PassiveVariables(nvariables, initial_condition, rhs!)

Append `nvariables` passive scalar variables to a [`SemidiscretizationImplicit`](@ref).
The passive variables are integrated by the time integrator but do not affect the
physical residual. The function `initial_condition(semi_base, t)` returns the initial
passive values, and `rhs!(dq, u_physical, du_physical, semi, t)` fills their time
derivatives.
"""
struct PassiveVariables{InitialCondition, RHS}
    nvariables::Int
    initial_condition::InitialCondition
    rhs!::RHS
end

function PassiveVariables(nvariables, initial_condition, rhs!)
    nvariables >= 0 ||
        throw(ArgumentError("Expected a non-negative number of passive variables."))
    return PassiveVariables{typeof(initial_condition), typeof(rhs!)}(Int(nvariables),
                                                                     initial_condition,
                                                                     rhs!)
end

@doc raw"""
    PassiveVariablesBoundaryFlux1D()

Append two passive scalar variables that store the time-integrated numerical boundary
fluxes at the negative and positive boundaries of a one-dimensional scalar problem. The
stored fluxes use the same sign convention as the parabolic flux in the semidiscrete
operator.
"""
struct PassiveVariablesBoundaryFlux1D end

@doc raw"""
    SemidiscretizationImplicit{Semidiscretization, TemporalOperator, PassiveVariables}

A semidiscretization wrapper that augments a spatial semidiscretization with a temporal
operator, so that the semi-discrete problem is not restricted to the explicit form
``\mathrm{d}\mathbf{u}/\mathrm{d}t = R(\mathbf{u}, t)``. For example,
[`TemporalOperatorConstitutive`](@ref) enforces a map between state and evolved
variables, while [`TemporalOperatorCapacity`](@ref) applies a nodal capacity function to
the spatial residual. Passive scalar variables can be appended to the ODE state for
diagnostics that are integrated by the time integrator but do not affect the physical
residual.
"""
struct SemidiscretizationImplicit{Semidiscretization <: Trixi.AbstractSemidiscretization,
                                  TemporalOperator <: AbstractTemporalOperator,
                                  PassiveVariables} <:
       Trixi.AbstractSemidiscretization
    semi_base::Semidiscretization
    operator_temporal::TemporalOperator
    passive_variables::PassiveVariables
end

function SemidiscretizationImplicit(semi_base::Trixi.AbstractSemidiscretization,
                                    operator_temporal::AbstractTemporalOperator)
    return SemidiscretizationImplicit(semi_base, operator_temporal,
                                      NoPassiveVariables())
end

function Base.show(io::IO, semi::SemidiscretizationImplicit)
    @nospecialize semi # reduce precompilation time
    print(io, "SemidiscretizationImplicit(")
    print(io, semi.semi_base)
    print(io, ", ", semi.operator_temporal |> typeof |> nameof)
    print(io, ", ", semi.passive_variables |> typeof |> nameof)
    print(io, ")")
    return nothing
end

@doc raw"""
    TemporalOperatorConstitutive(state_to_evolved; evolved_to_state = nothing)

Temporal operator for a [`SemidiscretizationImplicit`](@ref) that takes the form
```math
\frac{\mathrm{d}\mathbf{u}_\mathrm{evolved}}{\mathrm{d} t} =
R(\mathbf{u}_\mathrm{state}, t), \qquad
0 = \mathbf{u}_\mathrm{evolved} - S(\mathbf{u}_\mathrm{state}),
```
where ``R`` is the operator associated with the spatial discretization `semi_base`, and
``S`` is `state_to_evolved`, which maps the state variable(s) to the evolved variable(s).
The optional `evolved_to_state` inverse is required by AMR to reconstruct the algebraic
state after transferring the evolved block. In the case of the Richards equation, for
example, the state variable is the pressure head, and the evolved variable is the water
content, which are related by the water retention curve (i.e., [`water_content`](@ref)).
HydroTrixi.jl then uses SciML's DAE solvers to integrate the resulting system of
equations based on the following mass-matrix form:
```math
\begin{bmatrix} I & 0 \\ 0 & 0 \end{bmatrix} \frac{\mathrm{d}}{\mathrm{d} t}
\begin{bmatrix} \mathbf{u}_\mathrm{evolved} \\ \mathbf{u}_\mathrm{state} \end{bmatrix} =
\begin{bmatrix}
R(\mathbf{u}_\mathrm{state}, t) \\
\mathbf{u}_\mathrm{evolved} - S(\mathbf{u}_\mathrm{state})
\end{bmatrix}.
"""
struct TemporalOperatorConstitutive{StateToEvolved, EvolvedToState} <:
       AbstractTemporalOperator
    state_to_evolved::StateToEvolved
    evolved_to_state::EvolvedToState
end

function TemporalOperatorConstitutive(state_to_evolved; evolved_to_state = nothing)
    return TemporalOperatorConstitutive(state_to_evolved, evolved_to_state)
end

@doc raw"""
    TemporalOperatorCapacity(capacity_function; transfer_variables, transfer_to_state)

Temporal operator for a [`SemidiscretizationImplicit`](@ref) that stores the state
variable directly and applies a nodal capacity function to the spatial residual,
```math
\frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
= C(\mathbf{u})^{-1} R(\mathbf{u}, t).
```
For the pressure-head form of the Richards equation, the state variable is ``\psi`` and
``C(\psi) = \mathrm{d}\theta / \mathrm{d}\psi``. The capacity must be strictly positive
at all nodal states. The optional AMR transfer maps convert the state to the transferred
variable before mesh adaptation and reconstruct the state afterwards.
"""
struct TemporalOperatorCapacity{CapacityFunction, TransferVariables, TransferToState} <:
       AbstractTemporalOperator
    capacity_function::CapacityFunction
    transfer_variables::TransferVariables
    transfer_to_state::TransferToState
end

# Default AMR transfer keeps the stored state variable unchanged
@inline amr_transfer_identity(u, equations) = u

function TemporalOperatorCapacity(capacity_function;
                                  transfer_variables = amr_transfer_identity,
                                  transfer_to_state = amr_transfer_identity)
    return TemporalOperatorCapacity(capacity_function, transfer_variables,
                                    transfer_to_state)
end

print_temporal_operator_summary(io::IO, ::AbstractTemporalOperator) = nothing

function print_temporal_operator_summary(io::IO,
                                         operator_temporal::TemporalOperatorConstitutive)
    Trixi.summary_line(io, "state to evolved", operator_temporal.state_to_evolved)
    return Trixi.summary_line(io, "evolved to state",
                              operator_temporal.evolved_to_state)
end

function print_temporal_operator_summary(io::IO,
                                         operator_temporal::TemporalOperatorCapacity)
    Trixi.summary_line(io, "capacity function", operator_temporal.capacity_function)
    Trixi.summary_line(io, "transfer variables",
                       operator_temporal.transfer_variables)
    return Trixi.summary_line(io, "transfer to state",
                              operator_temporal.transfer_to_state)
end

@inline passive_variable_count(::NoPassiveVariables) = 0
@inline function passive_variable_count(passive_variables::PassiveVariables)
    return passive_variables.nvariables
end
@inline passive_variable_count(::PassiveVariablesBoundaryFlux1D) = 2

print_passive_variables_summary(io::IO, ::NoPassiveVariables) = nothing

function print_passive_variables_summary(io::IO, passive_variables)
    return Trixi.summary_line(io, "passive variables",
                              passive_variable_count(passive_variables))
end

# Wrapper to drive dispatch based on the temporal operator type on methods that take
# mesh, equations, solver, and cache as separate arguments.
struct CacheImplicit{Cache, TemporalOperator <: AbstractTemporalOperator, PassiveVariables}
    cache_base::Cache
    operator_temporal::TemporalOperator
    passive_variables::PassiveVariables
end

# Most properties are inherited from the base semidiscretization.
@inline function Base.getproperty(semi::SemidiscretizationImplicit, field::Symbol)
    if field === :performance_counter
        return getproperty(getfield(semi, :semi_base), field)
    end
    return getfield(semi, field)
end

@inline function Base.getproperty(cache::CacheImplicit, field::Symbol)
    if field === :cache_base || field === :operator_temporal ||
       field === :passive_variables
        return getfield(cache, field)
    end
    return getproperty(getfield(cache, :cache_base), field)
end

@inline Base.ndims(semi::SemidiscretizationImplicit) = ndims(semi.semi_base)
@inline Base.real(semi::SemidiscretizationImplicit) = real(semi.semi_base)

function Base.show(io::IO, ::MIME"text/plain", semi::SemidiscretizationImplicit)
    @nospecialize semi # reduce precompilation time

    if get(io, :compact, false)
        show(io, semi)
    else
        semi_base = semi.semi_base

        Trixi.summary_header(io, "SemidiscretizationImplicit")
        Trixi.summary_line(io, "#spatial dimensions", ndims(semi))
        Trixi.summary_line(io, "mesh", semi_base.mesh)
        Trixi.summary_line(io, "equations", semi_base.equations |> typeof |> nameof)
        Trixi.summary_line(io, "initial condition", semi_base.initial_condition)
        print_boundary_conditions_summary(io, semi_base.boundary_conditions)
        Trixi.summary_line(io, "source terms", semi_base.source_terms)
        Trixi.summary_line(io, "solver", semi_base.solver |> typeof |> nameof)
        Trixi.summary_line(io, "parabolic solver",
                           semi_base.solver_parabolic |> typeof |> nameof)
        Trixi.summary_line(io, "temporal operator",
                           semi.operator_temporal |> typeof |> nameof)
        print_temporal_operator_summary(io, semi.operator_temporal)
        print_passive_variables_summary(io, semi.passive_variables)
        Trixi.summary_line(io, "total #DOFs per field", Trixi.ndofsglobal(semi))
        Trixi.summary_footer(io)
    end
end

@inline function passive_variable_count(semi::SemidiscretizationImplicit)
    return passive_variable_count(semi.passive_variables)
end

@inline function passive_variable_count(cache::CacheImplicit)
    return passive_variable_count(cache.passive_variables)
end

# The full ODE state stores the physical DAE state first and passive scalars last
@inline function physical_variable_view(u_ode, semi::SemidiscretizationImplicit)
    n_passive = passive_variable_count(semi)
    return @view(u_ode[1:(length(u_ode) - n_passive)])
end

# Cache wrappers carry the passive layout needed by analysis integrations
@inline function physical_variable_view(u_ode, cache::CacheImplicit)
    n_passive = passive_variable_count(cache)
    return @view(u_ode[1:(length(u_ode) - n_passive)])
end

@doc raw"""
    passive_variable_view(u_ode, semi::SemidiscretizationImplicit)

Return a view of the passive scalar variables appended to the ODE state of `semi`.
"""
@inline function passive_variable_view(u_ode, semi::SemidiscretizationImplicit)
    n_passive = passive_variable_count(semi)
    n_passive == 0 && return @view(u_ode[1:0])
    return @view(u_ode[(length(u_ode) - n_passive + 1):length(u_ode)])
end

@doc raw"""
    passive_variables(u_ode, semi::SemidiscretizationImplicit)

Return a copy of the passive scalar variables appended to the ODE state of `semi`.
"""
function passive_variables(u_ode, semi::SemidiscretizationImplicit)
    return collect(passive_variable_view(u_ode, semi))
end

@doc raw"""
    boundary_flux_integrals(u_ode, semi::SemidiscretizationImplicit)

Return the integrated negative- and positive-boundary fluxes stored by
[`PassiveVariablesBoundaryFlux1D`](@ref) as a named tuple with fields `x_neg` and
`x_pos`.
"""
function boundary_flux_integrals(u_ode, semi::SemidiscretizationImplicit)
    semi.passive_variables isa PassiveVariablesBoundaryFlux1D ||
        throw(ArgumentError("Boundary flux integrals require " *
                            "`PassiveVariablesBoundaryFlux1D`."))
    passive_values = passive_variable_view(u_ode, semi)
    return (; x_neg = passive_values[1], x_pos = passive_values[2])
end

# Temporal operators split only the physical DAE state
@inline function evolved_variable_view(u_physical, ::TemporalOperatorCapacity)
    return u_physical
end

@inline function state_variable_view(u_physical, ::TemporalOperatorCapacity)
    return u_physical
end

@inline function evolved_variable_view(u_physical,
                                       ::TemporalOperatorConstitutive)
    return @view(u_physical[1:(length(u_physical) ÷ 2)])
end

@inline function state_variable_view(u_physical,
                                     ::TemporalOperatorConstitutive)
    return @view(u_physical[(length(u_physical) ÷ 2 + 1):end])
end

@inline function evolved_variable_view(u_ode, semi::SemidiscretizationImplicit)
    return evolved_variable_view(physical_variable_view(u_ode, semi),
                                 semi.operator_temporal)
end

@inline function state_variable_view(u_ode, semi::SemidiscretizationImplicit)
    return state_variable_view(physical_variable_view(u_ode, semi),
                               semi.operator_temporal)
end

# Error analysis compares state variables for both implicit forms
function Trixi.calc_error_norms(func, u_ode, t, analyzer,
                                semi::SemidiscretizationImplicit, cache_analysis)
    state_variable = state_variable_view(u_ode, semi)
    return Trixi.calc_error_norms(func, state_variable, t, analyzer, semi.semi_base,
                                  cache_analysis)
end

# Analysis callbacks operate on the physical state without passive variables
@inline function Trixi.wrap_array(u_ode::AbstractVector, mesh::Trixi.AbstractMesh,
                                  equations, dg::Trixi.DGSEM,
                                  cache::CacheImplicit{<:Any, <:TemporalOperatorCapacity})
    u_physical = physical_variable_view(u_ode, cache)
    return invoke(Trixi.wrap_array,
                  Tuple{AbstractVector, Trixi.AbstractMesh, Any, Trixi.DGSEM, Any},
                  u_physical, mesh, equations, dg, cache.cache_base)
end

# Constitutive analysis uses the evolved block as the conserved variable
@inline function Trixi.wrap_array(u_ode::AbstractVector, mesh::Trixi.AbstractMesh,
                                  equations, dg::Trixi.DGSEM,
                                  cache::CacheImplicit{<:Any,
                                                       <:TemporalOperatorConstitutive})
    evolved_variable = evolved_variable_view(physical_variable_view(u_ode, cache),
                                             cache.operator_temporal)
    return invoke(Trixi.wrap_array,
                  Tuple{AbstractVector, Trixi.AbstractMesh, Any, Trixi.DGSEM, Any},
                  evolved_variable, mesh, equations, dg, cache.cache_base)
end

# Default operator hooks correspond to the standard semidiscrete form `∂_t u = R(u, t)`.
@inline function nvariables_total(::AbstractTemporalOperator,
                                  semi_base::Trixi.AbstractSemidiscretization)
    return Trixi.nvariables(semi_base)
end

@inline check_ode_state(u_ode, ::AbstractTemporalOperator) = nothing

@inline function rhs_implicit!(du_ode, u_ode,
                               ::AbstractTemporalOperator,
                               semi_base::Trixi.AbstractSemidiscretization, t)
    return Trixi.default_rhs(semi_base)(du_ode, u_ode, semi_base, t)
end

@inline function mass_matrix(u_ode, ::AbstractTemporalOperator,
                             semi_base::Trixi.AbstractSemidiscretization)
    return Diagonal(ones(eltype(u_ode), length(u_ode)))
end

function rhs_implicit!(du_ode, u_ode, operator_temporal::TemporalOperatorCapacity,
                       semi_base::Trixi.AbstractSemidiscretization, t)
    Trixi.default_rhs(semi_base)(du_ode, u_ode, semi_base, t)
    (; equations) = semi_base
    capacity_function = operator_temporal.capacity_function

    # Apply the nodal capacity after assembling the spatial residual
    @inbounds for i in eachindex(du_ode, u_ode)
        du_ode[i] /= capacity_function(u_ode[i], equations)
    end

    return nothing
end

@inline function nvariables_total(::TemporalOperatorConstitutive,
                                  semi_base::Trixi.AbstractSemidiscretization)
    return 2 * Trixi.nvariables(semi_base)
end

@inline function check_ode_state(u_ode, ::TemporalOperatorConstitutive)
    iseven(length(u_ode)) ||
        throw(ArgumentError("Expected an even number of implicit degrees of freedom."))
    return nothing
end

function check_ode_state(u_ode, semi::SemidiscretizationImplicit)
    n_passive = passive_variable_count(semi)
    length(u_ode) >= n_passive ||
        throw(ArgumentError("Expected at least $n_passive passive variables."))
    return check_ode_state(physical_variable_view(u_ode, semi),
                           semi.operator_temporal)
end

function rhs_implicit!(du_ode, u_ode, operator_temporal::TemporalOperatorConstitutive,
                       semi_base::Trixi.AbstractSemidiscretization, t)
    evolved_variable = evolved_variable_view(u_ode, operator_temporal)
    state_variable = state_variable_view(u_ode, operator_temporal)
    evolved_variable_rhs = evolved_variable_view(du_ode, operator_temporal)
    state_variable_rhs = state_variable_view(du_ode, operator_temporal)

    # First block of du_ode: R(u_state, t)
    Trixi.default_rhs(semi_base)(evolved_variable_rhs, state_variable, semi_base, t)
    (; equations) = semi_base

    # Second block of du_ode: u_evolved - state_to_evolved(u_state)
    state_to_evolved = operator_temporal.state_to_evolved
    @inbounds for i in eachindex(state_variable_rhs, evolved_variable, state_variable)
        state_variable_rhs[i] = evolved_variable[i] -
                                state_to_evolved(state_variable[i], equations)
    end

    return nothing
end

function mass_matrix(u_ode, ::TemporalOperatorConstitutive,
                     semi_base::Trixi.AbstractSemidiscretization)
    half = length(u_ode) ÷ 2
    diagonal_entries = zeros(eltype(u_ode), length(u_ode))
    @inbounds diagonal_entries[1:half] .= one(eltype(u_ode))
    return Diagonal(diagonal_entries)
end

function passive_initial_values(::NoPassiveVariables, semi_base, t, RealT)
    return RealT[]
end

function passive_initial_values(passive_variables::PassiveVariables, semi_base, t,
                                RealT)
    values = passive_variables.initial_condition(semi_base, t)
    length(values) == passive_variable_count(passive_variables) ||
        throw(ArgumentError("Passive initial condition returned $(length(values)) " *
                            "values, expected " *
                            "$(passive_variable_count(passive_variables))."))

    passive_values = Vector{RealT}(undef, passive_variable_count(passive_variables))
    @inbounds for i in eachindex(passive_values)
        passive_values[i] = values[i]
    end

    return passive_values
end

function passive_initial_values(::PassiveVariablesBoundaryFlux1D, semi_base, t, RealT)
    return zeros(RealT, 2)
end

function set_passive_initial_values!(u_ode, semi::SemidiscretizationImplicit, t)
    passive_values = passive_initial_values(semi.passive_variables, semi.semi_base, t,
                                            eltype(u_ode))
    passive_variable_view(u_ode, semi) .= passive_values
    return nothing
end

function rhs_passive!(du_passive, u_physical, du_physical, ::NoPassiveVariables,
                      semi::SemidiscretizationImplicit, t)
    return nothing
end

function rhs_passive!(du_passive, u_physical, du_physical,
                      passive_variables::PassiveVariables,
                      semi::SemidiscretizationImplicit, t)
    return passive_variables.rhs!(du_passive, u_physical, du_physical, semi, t)
end

function rhs_passive!(du_passive, u_physical, du_physical,
                      ::PassiveVariablesBoundaryFlux1D,
                      semi::SemidiscretizationImplicit, t)
    flux_neg, flux_pos = boundary_fluxes_1d(semi.semi_base)
    du_passive[1] = flux_neg
    du_passive[2] = flux_pos
    return nothing
end

function boundary_fluxes_1d(semi_base::Trixi.AbstractSemidiscretization)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi_base)
    ndims(mesh) == 1 ||
        throw(ArgumentError("Boundary flux passive variables require a 1D mesh."))
    Trixi.nvariables(equations) == 1 ||
        throw(ArgumentError("Boundary flux passive variables require one variable."))

    n_boundaries_per_direction = cache.boundaries.n_boundaries_per_direction
    length(n_boundaries_per_direction) == 2 ||
        throw(ArgumentError("Expected two boundary directions in 1D."))
    n_boundaries_per_direction[1] == 1 && n_boundaries_per_direction[2] == 1 ||
        throw(ArgumentError("Boundary flux passive variables require non-periodic " *
                            "1D boundaries."))

    # Boundary fluxes are read from the cache filled by the physical RHS evaluation
    surface_flux_values = cache.elements.surface_flux_values
    lasts = accumulate(+, n_boundaries_per_direction)
    firsts = lasts - n_boundaries_per_direction .+ 1
    boundary_neg = firsts[1]
    boundary_pos = firsts[2]
    element_neg = cache.boundaries.neighbor_ids[boundary_neg]
    element_pos = cache.boundaries.neighbor_ids[boundary_pos]

    return surface_flux_values[1, 1, element_neg],
           surface_flux_values[1, 2, element_pos]
end

function mass_matrix(u_ode, semi::SemidiscretizationImplicit)
    u_physical = physical_variable_view(u_ode, semi)
    physical_mass_matrix = mass_matrix(u_physical, semi.operator_temporal,
                                       semi.semi_base)
    n_passive = passive_variable_count(semi)
    n_passive == 0 && return physical_mass_matrix

    # Passive scalar variables are differential variables
    passive_diagonal = ones(eltype(u_ode), n_passive)
    return Diagonal(vcat(physical_mass_matrix.diag, passive_diagonal))
end

@inline function Trixi.mesh_equations_solver_cache(semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return mesh, equations, solver,
           CacheImplicit(cache, semi.operator_temporal, semi.passive_variables)
end

@inline Trixi.default_rhs(::SemidiscretizationImplicit) = rhs_implicit!

# Number of variables depends on the temporal operator type, so we dispatch on that
@inline function Trixi.nvariables(semi::SemidiscretizationImplicit)
    return nvariables_total(semi.operator_temporal, semi.semi_base)
end

function implicit_physical_coefficients(t, semi_base::Trixi.AbstractSemidiscretization,
                                        ::TemporalOperatorCapacity)
    return Trixi.compute_coefficients(t, semi_base)
end

function implicit_physical_coefficients(t, semi_base::Trixi.AbstractSemidiscretization,
                                        operator_temporal::TemporalOperatorConstitutive)
    coefficients_state = Trixi.compute_coefficients(t, semi_base)
    coefficients_evolved = similar(coefficients_state)
    equations = semi_base.equations
    state_to_evolved = operator_temporal.state_to_evolved

    @inbounds for i in eachindex(coefficients_evolved, coefficients_state)
        coefficients_evolved[i] = state_to_evolved(coefficients_state[i], equations)
    end

    return vcat(coefficients_evolved, coefficients_state)
end

function implicit_physical_coefficients!(u_physical, t,
                                         semi_base::Trixi.AbstractSemidiscretization,
                                         ::TemporalOperatorCapacity)
    return Trixi.compute_coefficients!(u_physical, t, semi_base)
end

function implicit_physical_coefficients!(u_physical, t,
                                         semi_base::Trixi.AbstractSemidiscretization,
                                         operator_temporal::TemporalOperatorConstitutive)
    evolved_variable = evolved_variable_view(u_physical, operator_temporal)
    state_variable = state_variable_view(u_physical, operator_temporal)
    Trixi.compute_coefficients!(state_variable, t, semi_base)
    equations = semi_base.equations
    state_to_evolved = operator_temporal.state_to_evolved

    @inbounds for i in eachindex(evolved_variable, state_variable)
        evolved_variable[i] = state_to_evolved(state_variable[i], equations)
    end

    return nothing
end

function Trixi.compute_coefficients(t, semi::SemidiscretizationImplicit)
    coefficients_physical = implicit_physical_coefficients(t, semi.semi_base,
                                                           semi.operator_temporal)
    if passive_variable_count(semi) == 0
        return coefficients_physical
    end

    # Passive initial values are appended after the physical DAE state
    coefficients_passive = passive_initial_values(semi.passive_variables, semi.semi_base,
                                                  t,
                                                  eltype(coefficients_physical))
    coefficients_ode = vcat(coefficients_physical, coefficients_passive)
    record_mass_bias_initial_storage!(coefficients_ode, semi)
    return coefficients_ode
end

function Trixi.compute_coefficients!(u_ode, t, semi::SemidiscretizationImplicit)
    check_ode_state(u_ode, semi)

    u_physical = physical_variable_view(u_ode, semi)
    implicit_physical_coefficients!(u_physical, t, semi.semi_base,
                                    semi.operator_temporal)

    if passive_variable_count(semi) > 0
        # Passive initial values are only written when a tail block exists
        set_passive_initial_values!(u_ode, semi, t)
        record_mass_bias_initial_storage!(u_ode, semi)
    end
    return nothing
end

function rhs_implicit!(du_ode, u_ode, semi::SemidiscretizationImplicit, t)
    u_physical = physical_variable_view(u_ode, semi)
    du_physical = physical_variable_view(du_ode, semi)
    du_passive = passive_variable_view(du_ode, semi)

    # The physical residual is independent of the passive diagnostic variables
    rhs_implicit!(du_physical, u_physical, semi.operator_temporal, semi.semi_base, t)

    # Passive variables are filled after the physical RHS has updated solver caches
    rhs_passive!(du_passive, u_physical, du_physical, semi.passive_variables, semi, t)
    return nothing
end

function Trixi.semidiscretize(semi::SemidiscretizationImplicit, tspan;
                              reset_threads = true)
    if reset_threads
        Trixi.Polyester.reset_threads!()
    end

    u0_ode = Trixi.compute_coefficients(first(tspan), semi)
    check_ode_state(u0_ode, semi)

    mass_matrix_implicit = mass_matrix(u0_ode, semi)
    ode_function_type = SciMLBase.ODEFunction{true, SciMLBase.FullSpecialize}
    ode_function = ode_function_type(rhs_implicit!;
                                     mass_matrix = mass_matrix_implicit)
    return SciMLBase.ODEProblem{true, SciMLBase.FullSpecialize}(ode_function,
                                                                u0_ode,
                                                                tspan,
                                                                semi)
end
