@doc raw"""
    AbstractTemporalOperator

Abstract supertype for temporal formulations used by
[`SemidiscretizationImplicit`](@ref). A temporal operator determines the state
layout, residual construction, mass matrix, initial coefficients, and adaptive mesh
refinement reconstruction.
"""
abstract type AbstractTemporalOperator end

@doc raw"""
    AbstractPassiveVariables

Abstract supertype for passive variable configurations used for diagnostics by
[`SemidiscretizationImplicit`](@ref). Passive variables refer to additional scalar
variables appended to the physical ODE or DAE state that are integrated by the time
integrator but do not affect the physical residual. For example, a passive variable can
store the time-integrated numerical flux at a boundary for diagnostic purposes (see
[`PassiveVariablesBoundaryFlux1D`](@ref)).
"""
abstract type AbstractPassiveVariables end

@doc raw"""
    NoPassiveVariables()

Passive variable configuration for implicit semidiscretizations without appended
diagnostic variables.
"""
struct NoPassiveVariables <: AbstractPassiveVariables end

@doc raw"""
    PassiveVariablesBoundaryFlux1D()

Append two passive scalar variables that store the time-integrated numerical boundary
fluxes at the negative and positive boundaries of a one-dimensional scalar problem,
following the solver flux output method (SFOM) proposed by Ireson et al. (2023). For a
Richards column on ``[0,L]``, measuring positive downwards, the stored values are
``\int_{t_0}^t \hat{f}_0(\tau)\,\mathrm{d}\tau`` and
``\int_{t_0}^t \hat{f}_K(\tau)\,\mathrm{d}\tau``, where ``\hat{f}_0`` and ``\hat{f}_K`` are
the numerical fluxes at the top and bottom ends of the column.

# References
- Ireson, A. M., Spiteri, R. J., Clark, M. P., Mathias, S. A. (2023).
  A simple, efficient, mass-conservative approach to solving Richards'
  equation (openRE, v1.0). *Geoscientific Model Development*, 16, 659-677.
  [DOI: 10.5194/gmd-16-659-2023](https://doi.org/10.5194/gmd-16-659-2023)
"""
struct PassiveVariablesBoundaryFlux1D <: AbstractPassiveVariables end

@doc raw"""
    SemidiscretizationImplicit{Semidiscretization, TemporalOperator,
                               PassiveVariables}
    SemidiscretizationImplicit(semi_base, operator_temporal,
                               passive_variables = NoPassiveVariables())

A semidiscretization wrapper for the constant mass-matrix system
```math
\boldsymbol{A}\dot{\boldsymbol{y}}(t) = \boldsymbol{\mathcal{F}}(\boldsymbol{y}(t),t).
```
arising from a physical system that may not be expressed explicitly in terms of the time
derivative of the state variables. The `TemporalOperator` type may be
[`TemporalOperatorStandard`](@ref), [`TemporalOperatorConstitutive`](@ref), or
[`TemporalOperatorCapacity`](@ref), which determine the structure of the temporal mass
matrix and the partitioning of the state variables, and correspond to time-derivative
terms in PDEs of the form ``\partial_t u = R(u, t)``,
``\partial_t \vartheta(u) = R(u, t)``, or ``C(u)\partial_t u = R(u, t)``, respectively,
where $R$ denotes a generic spatial operator. The `PassiveVariables` type may be
[`NoPassiveVariables`](@ref) or [`PassiveVariablesBoundaryFlux1D`](@ref), the latter of
which appends two passive scalar variables to the ODE state that store the time-integrated
boundary fluxes for a one-dimensional scalar problem.

!!! note
    The constant temporal mass matrix ``\boldsymbol{A}`` is distinct from the spatial
    discretization mass matrix, which is handled by `semi_base` inside
    ``\boldsymbol{\mathcal{R}}``.
"""
struct SemidiscretizationImplicit{Semidiscretization <: Trixi.AbstractSemidiscretization,
                                  TemporalOperator <: AbstractTemporalOperator,
                                  PassiveVariables <: AbstractPassiveVariables} <:
       Trixi.AbstractSemidiscretization
    semi_base::Semidiscretization
    operator_temporal::TemporalOperator
    passive_variables::PassiveVariables
end

function SemidiscretizationImplicit(semi_base::Trixi.AbstractSemidiscretization,
                                    operator_temporal::AbstractTemporalOperator)
    return SemidiscretizationImplicit(semi_base, operator_temporal, NoPassiveVariables())
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
    TemporalOperatorStandard()

Temporal operator for the standard semidiscrete form
```math
\dot{\boldsymbol{u}}(t) = \boldsymbol{\mathcal{R}}(\boldsymbol{u}(t),t).
```
Before passive variables are appended, the [`SemidiscretizationImplicit`](@ref) state
satisfies ``\boldsymbol{y} = \boldsymbol{u}``,
``\boldsymbol{\mathcal{F}} = \boldsymbol{\mathcal{R}}``, and the mass matrix is the
identity.
"""
struct TemporalOperatorStandard <: AbstractTemporalOperator end

@doc raw"""
    TemporalOperatorConstitutive(state_to_evolved; evolved_to_state = nothing)

Temporal operator for a [`SemidiscretizationImplicit`](@ref) that takes the form
```math
\begin{bmatrix} I & 0 \\ 0 & 0 \end{bmatrix}
\frac{\mathrm{d}}{\mathrm{d}t}
\begin{bmatrix}
\boldsymbol{u}_\mathrm{evolved} \\ \boldsymbol{u}_\mathrm{state}
\end{bmatrix}
=
\begin{bmatrix}
\boldsymbol{\mathcal{R}}(\boldsymbol{u}_\mathrm{state},t) \\
\boldsymbol{u}_\mathrm{evolved} -
\boldsymbol{\vartheta}(\boldsymbol{u}_\mathrm{state})
\end{bmatrix}.
```
Here, ``\boldsymbol{\mathcal{R}}`` is the spatial operator, and
``\boldsymbol{\vartheta}`` is `state_to_evolved`, the generic constitutive map. Thus,
``\boldsymbol{y}`` contains distinct blocks ordered as evolved variables followed by
state variables. Passive variables, when present, are appended after both blocks. The
optional `evolved_to_state` inverse is required by adaptive mesh refinement to
reconstruct the state variables after transferring the evolved block.
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
variable directly and applies a nodal capacity function to the spatial operator,
```math
\dot{\boldsymbol{u}}(t) =
\boldsymbol{C}(\boldsymbol{u}(t))^{-1}
\boldsymbol{\mathcal{R}}(\boldsymbol{u}(t),t).
```
The capacity must be strictly positive at all nodal states. Before passive variables are
appended, the stored state satisfies ``\boldsymbol{y} = \boldsymbol{u}``, the residual is
``\boldsymbol{\mathcal{F}} = \boldsymbol{C}^{-1}\boldsymbol{\mathcal{R}}``, and the mass
matrix is the identity. The optional adaptive mesh refinement transfer maps convert the
state variables to the transferred variables before mesh adaptation and reconstruct the
state afterwards.
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
    return Trixi.summary_line(io, "evolved to state", operator_temporal.evolved_to_state)
end

function print_temporal_operator_summary(io::IO,
                                         operator_temporal::TemporalOperatorCapacity)
    Trixi.summary_line(io, "capacity function", operator_temporal.capacity_function)
    Trixi.summary_line(io, "transfer variables", operator_temporal.transfer_variables)
    return Trixi.summary_line(io, "transfer to state", operator_temporal.transfer_to_state)
end

@inline passive_variable_count(::NoPassiveVariables) = 0
@inline passive_variable_count(::PassiveVariablesBoundaryFlux1D) = 2

print_passive_variables_summary(io::IO, ::NoPassiveVariables) = nothing

function print_passive_variables_summary(io::IO, passive_variables)
    return Trixi.summary_line(io, "passive variables",
                              passive_variable_count(passive_variables))
end

# Wrapper to drive dispatch based on the temporal operator type on methods that take
# mesh, equations, solver, and cache as separate arguments.
struct CacheImplicit{Cache, TemporalOperator <: AbstractTemporalOperator,
                     PassiveVariables <: AbstractPassiveVariables}
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
    if field === :cache_base || field === :operator_temporal || field === :passive_variables
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

# Return the number of passive variables appended to the ODE state of `semi`. The total
# number of degrees of freedom is the sum of the physical and passive degrees of freedom.
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

# Return a view of passive diagnostic variables appended to the ODE state
@inline function passive_variable_view(u_ode, semi::SemidiscretizationImplicit)
    n_passive = passive_variable_count(semi)
    n_passive == 0 && return @view(u_ode[1:0])
    return @view(u_ode[(length(u_ode) - n_passive + 1):length(u_ode)])
end

# Return a copy of passive diagnostic variables appended to the ODE state
function passive_variables(u_ode, semi::SemidiscretizationImplicit)
    return collect(passive_variable_view(u_ode, semi))
end

# Return integrated negative- and positive-boundary fluxes stored as passive variables
function boundary_flux_integrals(u_ode, semi::SemidiscretizationImplicit)
    passive_values = passive_variable_view(u_ode, semi)
    return (; x_neg = passive_values[1], x_pos = passive_values[2])
end

# Standard and capacity operators use one stored vector for both variable roles
@inline function evolved_variable_view(u_physical,
                                       ::Union{TemporalOperatorStandard,
                                               TemporalOperatorCapacity})
    return u_physical
end

# Standard and capacity operators supply the stored vector to the spatial operator
@inline function state_variable_view(u_physical,
                                     ::Union{TemporalOperatorStandard,
                                             TemporalOperatorCapacity})
    return u_physical
end

@inline function evolved_variable_view(u_physical, ::TemporalOperatorConstitutive)
    return @view(u_physical[1:(length(u_physical) ÷ 2)])
end

@inline function state_variable_view(u_physical, ::TemporalOperatorConstitutive)
    return @view(u_physical[(length(u_physical) ÷ 2 + 1):end])
end

@inline function evolved_variable_view(u_ode, semi::SemidiscretizationImplicit)
    return evolved_variable_view(physical_variable_view(u_ode, semi),
                                 semi.operator_temporal)
end

@inline function state_variable_view(u_ode, semi::SemidiscretizationImplicit)
    return state_variable_view(physical_variable_view(u_ode, semi), semi.operator_temporal)
end

# Error analysis compares u for all temporal formulations
function Trixi.calc_error_norms(func, u_ode, t, analyzer,
                                semi::SemidiscretizationImplicit, cache_analysis)
    state_variable = state_variable_view(u_ode, semi)
    return Trixi.calc_error_norms(func, state_variable, t, analyzer, semi.semi_base,
                                  cache_analysis)
end

# Standard and capacity analysis use u directly
@inline function Trixi.wrap_array(u_ode::AbstractVector, mesh::Trixi.AbstractMesh,
                                  equations, dg::Trixi.DGSEM,
                                  cache::CacheImplicit{<:Any,
                                                       <:Union{TemporalOperatorStandard,
                                                               TemporalOperatorCapacity}})
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
                                  semi_base)
    return Trixi.nvariables(semi_base)
end

@inline function rhs_implicit!(du_ode, u_ode, ::AbstractTemporalOperator,
                               semi_base, t)
    return Trixi.default_rhs(semi_base)(du_ode, u_ode, semi_base, t)
end

@inline function mass_matrix(u_ode, ::AbstractTemporalOperator,
                             semi_base)
    return Diagonal(ones(eltype(u_ode), length(u_ode)))
end

function rhs_implicit!(du_ode, u_ode, operator_temporal::TemporalOperatorCapacity,
                       semi_base, t)
    Trixi.default_rhs(semi_base)(du_ode, u_ode, semi_base, t)
    (; equations) = semi_base
    capacity_function = operator_temporal.capacity_function

    # Apply the nodal capacity after assembling the spatial residual
    @inbounds for i in eachindex(du_ode, u_ode)
        du_ode[i] /= capacity_function(u_ode[i], equations)
    end

    return nothing
end

# The constitutive operator stores evolved variables first and state variables second
@inline function nvariables_total(::TemporalOperatorConstitutive,
                                  semi_base)
    return 2 * Trixi.nvariables(semi_base)
end

function rhs_implicit!(du_ode, u_ode, operator_temporal::TemporalOperatorConstitutive,
                       semi_base, t)
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
                     semi_base)
    half = length(u_ode) ÷ 2
    diagonal_entries = zeros(eltype(u_ode), length(u_ode))
    @inbounds diagonal_entries[1:half] .= one(eltype(u_ode))
    return Diagonal(diagonal_entries)
end

function passive_initial_values(::NoPassiveVariables, semi_base, t, RealT)
    return RealT[]
end

function passive_initial_values(::PassiveVariablesBoundaryFlux1D, semi_base, t, RealT)
    return zeros(RealT, 2)
end

function rhs_passive!(du_passive, u_physical, du_physical, ::NoPassiveVariables,
                      semi, t)
    return nothing
end

function rhs_passive!(du_passive, u_physical, du_physical, ::PassiveVariablesBoundaryFlux1D,
                      semi::SemidiscretizationImplicit{<:Trixi.SemidiscretizationParabolic{<:Trixi.AbstractMesh{1},
                                                                                           <:Trixi.AbstractEquationsParabolic{1,
                                                                                                                              1},
                                                                                           <:Any,
                                                                                           <:NamedTuple{(:x_neg,
                                                                                                         :x_pos)}}},
                      t)
    cache = semi.semi_base.cache
    n_boundaries_per_direction = cache.boundaries.n_boundaries_per_direction

    # Boundary fluxes are read from the cache filled by the physical RHS evaluation
    surface_flux_values = cache.elements.surface_flux_values
    lasts = accumulate(+, n_boundaries_per_direction)
    firsts = lasts - n_boundaries_per_direction .+ 1
    boundary_neg = firsts[1]
    boundary_pos = firsts[2]
    element_neg = cache.boundaries.neighbor_ids[boundary_neg]
    element_pos = cache.boundaries.neighbor_ids[boundary_pos]

    du_passive[1] = surface_flux_values[1, 1, element_neg]
    du_passive[2] = surface_flux_values[1, 2, element_pos]
    return nothing
end

function mass_matrix(u_ode, semi::SemidiscretizationImplicit)
    u_physical = physical_variable_view(u_ode, semi)
    physical_mass_matrix = mass_matrix(u_physical, semi.operator_temporal, semi.semi_base)
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

# Standard and capacity operators initialize u directly
function implicit_physical_coefficients(t, semi_base,
                                        ::Union{TemporalOperatorStandard,
                                                TemporalOperatorCapacity})
    return Trixi.compute_coefficients(t, semi_base)
end

function implicit_physical_coefficients(t, semi_base,
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

# Standard and capacity operators write initial data directly to their shared vector
function implicit_physical_coefficients!(u_physical, t, semi_base,
                                         ::Union{TemporalOperatorStandard,
                                                 TemporalOperatorCapacity})
    return Trixi.compute_coefficients!(u_physical, t, semi_base)
end

function implicit_physical_coefficients!(u_physical, t,
                                         semi_base,
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
    coefficients_passive = passive_initial_values(semi.passive_variables, semi.semi_base, t,
                                                  eltype(coefficients_physical))
    coefficients_ode = vcat(coefficients_physical, coefficients_passive)
    record_mass_bias_initial_storage!(coefficients_ode, semi, semi.operator_temporal)
    return coefficients_ode
end

function Trixi.compute_coefficients!(u_ode, t, semi::SemidiscretizationImplicit)
    u_physical = physical_variable_view(u_ode, semi)
    implicit_physical_coefficients!(u_physical, t, semi.semi_base, semi.operator_temporal)

    if passive_variable_count(semi) > 0
        # Passive initial values are only written when a tail block exists
        passive_values = passive_initial_values(semi.passive_variables, semi.semi_base,
                                                t, eltype(u_ode))
        passive_variable_view(u_ode, semi) .= passive_values
        record_mass_bias_initial_storage!(u_ode, semi, semi.operator_temporal)
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

@doc raw"""
    semidiscretize(semi::SemidiscretizationImplicit, tspan;
                   reset_threads = true, jacobian = SparseJacobian())

Construct a `SciMLBase.ODEProblem` for the constant mass-matrix system represented by a
semidiscretization of type [`SemidiscretizationImplicit`](@ref). The `jacobian` strategy
controls which Jacobian information HydroTrixi.jl supplies to SciML. See
[`DenseJacobian`](@ref) and [`SparseJacobian`](@ref) for the available strategies.

!!! warning
    `SparseJacobian()` supports only serial, nonperiodic, scalar, one-dimensional
    `TreeMesh` problems using a Lobatto-Legendre `DGSEM` and
    `ParabolicFormulationLocalDG` with any penalty parameter. The supported temporal
    operators are
    [`TemporalOperatorStandard`](@ref), [`TemporalOperatorCapacity`](@ref), and
    [`TemporalOperatorConstitutive`](@ref); the supported passive variable types are
    [`NoPassiveVariables`](@ref) and [`PassiveVariablesBoundaryFlux1D`](@ref). MPI
    execution and periodic meshes are rejected explicitly; other configurations are
    unsupported and fail through ordinary Julia dispatch.
"""
function Trixi.semidiscretize(semi::SemidiscretizationImplicit, tspan; reset_threads = true,
                              jacobian = SparseJacobian())
    if reset_threads
        Trixi.Polyester.reset_threads!
    end

    u0_ode = Trixi.compute_coefficients(first(tspan), semi)

    mass_matrix_implicit = mass_matrix(u0_ode, semi)
    jacobian_options_implicit = jacobian_options(jacobian, u0_ode, semi)
    rhs_implicit_cached = RHSImplicitCache()
    ode_function_type = SciMLBase.ODEFunction{true, SciMLBase.FullSpecialize}
    ode_function = ode_function_type(rhs_implicit_cached;
                                     mass_matrix = mass_matrix_implicit,
                                     jacobian_options_implicit...)
    return SciMLBase.ODEProblem{true, SciMLBase.FullSpecialize}(ode_function, u0_ode, tspan,
                                                                semi)
end
