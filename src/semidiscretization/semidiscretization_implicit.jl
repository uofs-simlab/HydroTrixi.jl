@doc raw"""
    AbstractTemporalOperator

Abstract supertype for temporal formulations used by [`SemidiscretizationImplicit`](@ref).
A temporal operator determines the state layout, residual construction, DAE mass matrix,
initial coefficients, and adaptive mesh refinement reconstruction.
"""
abstract type AbstractTemporalOperator end

@doc raw"""
    AbstractPassiveVariables

Abstract supertype for passive variable configurations used for diagnostics by
[`SemidiscretizationImplicit`](@ref). Passive variables refer to additional scalar variables
appended to the physical ODE or DAE state that are integrated by the time integrator but do
not affect the physical residual. For example, a passive variable can store a cumulative
numerical boundary flux for diagnostic purposes (see
[`PassiveVariablesBoundaryFlux1D`](@ref)).
"""
abstract type AbstractPassiveVariables end

@doc raw"""
    NoPassiveVariables()

Passive variable configuration for implicit semidiscretizations without appended diagnostic
variables.
"""
struct NoPassiveVariables <: AbstractPassiveVariables end

@doc raw"""
    PassiveVariablesBoundaryFlux1D()

Append two passive scalar variables that store the cumulative numerical boundary fluxes at
the negative and positive boundaries of a one-dimensional scalar problem, following the
solver flux output method (SFOM) proposed by Ireson et al. (2023). For a Richards column on
``[0,L]``, measuring depth positive downward, denote these variables by ``F_{\mathrm{T}}``
and ``F_{\mathrm{B}}``, satisfying

```math
\dot{F}_{\mathrm{T}}(t)=f_{\mathrm{T}}^\star(t), \quad
\dot{F}_{\mathrm{B}}(t)=f_{\mathrm{B}}^\star(t), \quad
F_{\mathrm{T}}(t^0)=F_{\mathrm{B}}(t^0)=0,
```

where ``f_{\mathrm{T}}^\star`` and ``f_{\mathrm{B}}^\star`` are the numerical fluxes at
the soil surface and bottom of the column, respectively. The returned named-tuple keys
`x_neg` and `x_pos` represent these top and bottom quantities in the Julia API.

# References

- Ireson, A. M., Spiteri, R. J., Clark, M. P., Mathias, S. A. (2023).
  A simple, efficient, mass-conservative approach to solving
  Richards' equation (openRE, v1.0). *Geoscientific Model Development*, 16, 659-677.
  [DOI: 10.5194/gmd-16-659-2023](https://doi.org/10.5194/gmd-16-659-2023)
"""
struct PassiveVariablesBoundaryFlux1D <: AbstractPassiveVariables end

@doc raw"""
    SemidiscretizationImplicit{Semidiscretization, TemporalOperator,
                               PassiveVariables, CacheParabolic}
    SemidiscretizationImplicit(semi_base, operator_temporal,
                               passive_variables = NoPassiveVariables())

A semidiscretization wrapper for the constant mass-matrix system
```math
\boldsymbol{M}\dot{\boldsymbol{y}}(t) =
\boldsymbol{\mathcal{F}}(\boldsymbol{y}(t),t).
```
arising from a physical system that may not be expressed explicitly in terms of the time
derivative of the state variables. The `TemporalOperator` type may be
[`TemporalOperatorStandard`](@ref), [`TemporalOperatorConstitutive`](@ref), or
[`TemporalOperatorCapacity`](@ref), which determine the structure of the DAE mass matrix and
the partitioning of the state variables, and correspond to time-derivative terms in PDEs of
the form ``\partial_t u = \mathcal{R}(u, t)``,
``\partial_t \vartheta(u) = \mathcal{R}(u, t)``, or
``c(u)\partial_t u = \mathcal{R}(u, t)``, respectively, where ``\mathcal{R}`` denotes a
generic spatial operator. The `PassiveVariables` type may be [`NoPassiveVariables`](@ref) or
[`PassiveVariablesBoundaryFlux1D`](@ref), the latter of which appends two passive scalar
variables to the ODE state that store the time-integrated boundary fluxes for a
one-dimensional scalar problem.

For one-dimensional parabolic `TreeMesh` discretizations with a Lobatto-Legendre `DGSEM`,
the wrapper adds boundary-flux storage without replacing `semi_base`. Its spatial operator
retains both the interior solution and flux at each boundary. Other base discretizations,
including hyperbolic and multidimensional parabolic systems, do not create a
`CacheParabolic1D`, and the implicit wrapper's `cache_parabolic` field is `nothing`. Their
spatial operator is evaluated by `Trixi.default_rhs(semi_base)`, using `semi_base.cache`
and, for parabolic systems, `semi_base.cache_parabolic`.

!!! note
    The constant DAE mass matrix ``\boldsymbol{M}``, constructed by `dae_mass_matrix`,
    multiplies the time derivative of the ODE or DAE state. The spatial discretization's
    mass matrix is handled separately by `semi_base` inside ``\boldsymbol{\mathcal{R}}``.
"""
struct SemidiscretizationImplicit{Semidiscretization <: Trixi.AbstractSemidiscretization,
                                  TemporalOperator <: AbstractTemporalOperator,
                                  PassiveVariables <: AbstractPassiveVariables,
                                  CacheParabolic} <:
       Trixi.AbstractSemidiscretization
    semi_base        ::Semidiscretization
    operator_temporal::TemporalOperator
    passive_variables::PassiveVariables
    cache_parabolic  ::CacheParabolic
end

function SemidiscretizationImplicit(semi_base::Trixi.AbstractSemidiscretization,
                                    operator_temporal::AbstractTemporalOperator,
                                    passive_variables::AbstractPassiveVariables)
    cache_parabolic = create_cache_parabolic_implicit(semi_base)
    return SemidiscretizationImplicit(semi_base, operator_temporal, passive_variables,
                                      cache_parabolic)
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

Temporal operator for the standard semi-discrete form
```math
\dot{\boldsymbol{u}}(t) = \boldsymbol{\mathcal{R}}(\boldsymbol{u}(t),t).
```
Before passive variables are appended, the [`SemidiscretizationImplicit`](@ref) state
satisfies ``\boldsymbol{y} = \boldsymbol{u}``,
``\boldsymbol{\mathcal{F}} = \boldsymbol{\mathcal{R}}``, and the DAE mass matrix is the
identity.
"""
struct TemporalOperatorStandard <: AbstractTemporalOperator end

@doc raw"""
    TemporalOperatorConstitutive(state_to_evolved; evolved_to_state = nothing,
                                 transfer_state = false)

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
Here, ``\boldsymbol{\mathcal{R}}`` is the spatial operator, and ``\boldsymbol{\vartheta}``
is `state_to_evolved`, the generic constitutive map. Thus, ``\boldsymbol{y}`` contains
distinct blocks ordered as evolved variables followed by state variables. Passive variables,
when present, are appended after both blocks. By default, adaptive mesh refinement transfers
the evolved block and reconstructs the state block with
`evolved_to_state(value, equations)`. Using AMR with `transfer_state = false` therefore
requires a callable `evolved_to_state`; the default `nothing` does not support this
reconstruction. With `transfer_state = true`, AMR transfers the state block and reconstructs
the evolved block with `state_to_evolved(value, equations)`.

For the mixed Richards formulation, these generic blocks are
``\boldsymbol{u}_{\mathrm{evolved}}=\boldsymbol{\Theta}`` and
``\boldsymbol{u}_{\mathrm{state}}=\boldsymbol{\Psi}``, while the spatial and constitutive
maps are ``\boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t)`` and
``\boldsymbol{\vartheta}(\boldsymbol{\Psi})``, respectively.
"""
struct TemporalOperatorConstitutive{StateToEvolved, EvolvedToState} <:
       AbstractTemporalOperator
    state_to_evolved::StateToEvolved
    evolved_to_state::EvolvedToState
    transfer_state::Bool
end

function TemporalOperatorConstitutive(state_to_evolved; evolved_to_state = nothing,
                                      transfer_state = false)
    return TemporalOperatorConstitutive(state_to_evolved, evolved_to_state,
                                        transfer_state)
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
``\boldsymbol{\mathcal{F}} = \boldsymbol{C}^{-1}\boldsymbol{\mathcal{R}}``, and the DAE mass
matrix is the identity. The optional adaptive mesh refinement transfer maps convert the
state variables to the transferred variables before mesh adaptation and reconstruct the
state afterwards.

For the pressure-head Richards formulation, ``\boldsymbol{u}=\boldsymbol{\Psi}`` and the
right-hand side is
``\boldsymbol{C}(\boldsymbol{\Psi})^{-1} \boldsymbol{\mathcal{R}}(\boldsymbol{\Psi},t)``.
"""
struct TemporalOperatorCapacity{CapacityFunction, TransferVariables, TransferToState} <:
       AbstractTemporalOperator
    capacity_function::CapacityFunction
    transfer_variables::TransferVariables
    transfer_to_state::TransferToState
end

# Default AMR transfer keeps the stored state variable unchanged
function TemporalOperatorCapacity(capacity_function;
                                  transfer_variables = (u, equations) -> u,
                                  transfer_to_state = (u, equations) -> u)
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

# Return the number of passive variables appended to the physical state. The total
# number of degrees of freedom is the sum of the physical and passive degrees of freedom.
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

# Read semi.performance_counter from semi.semi_base.performance_counter.
# Read all other properties directly from semi's own fields.
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

# The full ODE state stores the physical DAE state first and passive scalars last
@inline function physical_variable_view(u_ode,
                                        passive_variables::AbstractPassiveVariables)
    n_passive = passive_variable_count(passive_variables)
    return @view(u_ode[1:(length(u_ode) - n_passive)])
end

# Return a view of passive diagnostic variables appended to the ODE state
@inline function passive_variable_view(u_ode, semi::SemidiscretizationImplicit)
    n_passive = passive_variable_count(semi.passive_variables)
    # With no passive variables, the tail range is empty: (length(u_ode) + 1):length(u_ode)
    return @view(u_ode[(length(u_ode) - n_passive + 1):length(u_ode)])
end

# Return integrated negative- and positive-boundary fluxes stored as passive variables
function boundary_flux_integrals(u_ode, semi::SemidiscretizationImplicit)
    passive_values = passive_variable_view(u_ode, semi)
    return (; x_neg = passive_values[1], x_pos = passive_values[2])
end

@doc raw"""
    evolved_variable_block(u_physical, operator)
    evolved_variable_block(u_ode, semi::SemidiscretizationImplicit)

Return the block containing the evolved variables. For [`TemporalOperatorStandard`](@ref)
and [`TemporalOperatorCapacity`](@ref), this is the complete physical state. For
[`TemporalOperatorConstitutive`](@ref), this is the first half of the physical state.

The returned block aliases its input and must not be assumed to be an independent copy.
The `semi` method also excludes appended passive diagnostic variables.
"""
@inline function evolved_variable_block(u_physical,
                                        ::Union{TemporalOperatorStandard,
                                                TemporalOperatorCapacity})
    return u_physical
end

@doc raw"""
    state_variable_block(u_physical, operator)
    state_variable_block(u_ode, semi::SemidiscretizationImplicit)

Return the block containing the state variables supplied to the spatial operator. For
[`TemporalOperatorStandard`](@ref) and [`TemporalOperatorCapacity`](@ref), this is the
complete physical state. For [`TemporalOperatorConstitutive`](@ref), this is the second half
of the physical state.

The returned block aliases its input and must not be assumed to be an independent copy.
The `semi` method also excludes appended passive diagnostic variables.
"""
@inline function state_variable_block(u_physical,
                                      ::Union{TemporalOperatorStandard,
                                              TemporalOperatorCapacity})
    return u_physical
end

@inline function evolved_variable_block(u_physical, ::TemporalOperatorConstitutive)
    return @view(u_physical[1:(length(u_physical) ÷ 2)])
end

@inline function state_variable_block(u_physical, ::TemporalOperatorConstitutive)
    return @view(u_physical[(length(u_physical) ÷ 2 + 1):end])
end

@inline function evolved_variable_block(u_ode, semi::SemidiscretizationImplicit)
    u_physical = physical_variable_view(u_ode, semi.passive_variables)
    return evolved_variable_block(u_physical, semi.operator_temporal)
end

@inline function state_variable_block(u_ode, semi::SemidiscretizationImplicit)
    u_physical = physical_variable_view(u_ode, semi.passive_variables)
    return state_variable_block(u_physical, semi.operator_temporal)
end

# Evaluate error norms using state_variable_block, which contains pressure head for
# Richards in both temporal formulations. The function func selects the quantity
# compared with the reference solution.
function Trixi.calc_error_norms(func, u_ode, t, analyzer,
                                semi::SemidiscretizationImplicit, cache_analysis)
    state_variable = state_variable_block(u_ode, semi)
    semi_base = semi.semi_base
    mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi_base)
    GC.@preserve u_ode begin
        u = unsafe_wrap_block(state_variable, mesh, equations, dg, cache)
        return Trixi.calc_error_norms(func, u, t, analyzer, mesh, equations,
                                      semi_base.initial_condition, dg, cache,
                                      cache_analysis)
    end
end

"""
    unsafe_vector_alias(u_block)

Return a `Vector` sharing the entries and element type of `u_block`, without copying. An
input `Vector` is returned unchanged; a view of consecutive array entries is exposed through
a new vector alias.

For a view `u_block` into the ODE vector `u_ode`, the returned `u_vector` points into
`u_ode`'s memory but holds no Julia reference to `u_ode`. The garbage collector may free
this memory after the last direct use of `u_ode`, even while `u_vector` is still being
used. Preserve `u_ode` with `GC.@preserve` throughout this call and all uses of `u_vector`.
Do not resize `u_ode` during that time.

```julia
GC.@preserve u_ode begin
    u_block = state_variable_block(u_ode, semi)
    u_vector = unsafe_vector_alias(u_block)
    # All reading and writing through u_vector goes here.
end
```
"""
@inline function unsafe_vector_alias(u_block::Vector)
    return u_block
end

@inline function unsafe_vector_alias(u_block::SubArray{<:Any, 1, <:Array,
                                                      <:Tuple{<:AbstractUnitRange}, true})
    # Trixi's plain-array wrapping path needs concrete Vector storage.
    # The pointer includes u_block's offset; own = false prevents the returned Vector
    # from freeing the memory belonging to parent(u_block).
    return unsafe_wrap(Vector{eltype(u_block)}, pointer(u_block), length(u_block);
                       own = false)
end

"""
    unsafe_wrap_block(u_block, mesh, equations, dg, cache)

Wrap a `Vector` or a view of consecutive array entries in Trixi.jl's
`(variable, node..., element)` layout without copying, and without changing the element
type. Select the required block (e.g., `evolved_variable_block(u_ode, semi)` or
`state_variable_block(u_ode, semi)`) and exclude passive variables before calling this
function.

For a block selected from the ODE vector `u_ode`, the returned nodal array `u` does not
automatically keep `u_ode` alive. Therefore, you must explicitly keep `u_ode` alive
throughout the use of `u` to prevent the garbage collector from freeing its memory
prematurely. This also applies when `u_block` is `u_ode` itself. Preserve `u_ode` with
`GC.@preserve` throughout this call and all uses of `u` and do not resize `u_ode` during
that time, for example:

```julia
mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
GC.@preserve u_ode begin
    u_block = state_variable_block(u_ode, semi)
    u = unsafe_wrap_block(u_block, mesh, equations, dg, cache)
    # All reading and writing through the returned nodal array u goes here.
end
```
"""
@inline function unsafe_wrap_block(u_block::AbstractVector, mesh, equations, dg, cache)
    # Expose the block as concrete Vector storage for Trixi
    u_vector = unsafe_vector_alias(u_block)

    # Trixi.wrap_array returns a multidimensional PtrArray over u_vector's entries
    # when supported. If eltype(u_block) <: ForwardDiff.Dual, it returns a plain Array
    # instead, since LoopVectorization does not support dual numbers. Neither path copies
    # entries.
    return Trixi.wrap_array(u_vector, mesh, equations, dg, cache)
end

# Standard and capacity formulations expose the physical state for analysis, excluding
# passive variables. For Richards in capacity form, this state is pressure head.
# Use GC.@preserve on u_ode throughout this call and all uses of the returned array
# to keep u_ode's storage alive.
@inline function Trixi.wrap_array(u_ode::AbstractVector, mesh::Trixi.AbstractMesh,
                                  equations, dg::Trixi.DGSEM,
                                  cache::CacheImplicit{<:Any,
                                                       <:Union{TemporalOperatorStandard,
                                                               TemporalOperatorCapacity}})
    # Extract a view of the physical variables from the ODE vector (i.e., no passive vars)
    u_physical = physical_variable_view(u_ode, cache.passive_variables)

    # Make the physical-variable view into a vector without copying values, and then pass
    # the flat vector to Trixi.jl's Trixi.wrap_array, which converts it into the
    # multidimensional layout expected by the spatial operators.
    return unsafe_wrap_block(u_physical, mesh, equations, dg, cache.cache_base)
end

# Analysis with TemporalOperatorConstitutive uses the evolved block as the conserved
# variable. Use GC.@preserve on u_ode throughout this call and all uses of the returned
# array to prevent the garbage collector from freeing its memory prematurely.
@inline function Trixi.wrap_array(u_ode::AbstractVector, mesh::Trixi.AbstractMesh,
                                  equations, dg::Trixi.DGSEM,
                                  cache::CacheImplicit{<:Any,
                                                       <:TemporalOperatorConstitutive})
    # Extract a view of the physical variables from the ODE vector (i.e., no passive vars)
    u_physical = physical_variable_view(u_ode, cache.passive_variables)

    # Get a view of only the evolved (i.e., differential) variable block
    evolved_variable = evolved_variable_block(u_physical, cache.operator_temporal)

    # Make the evolved-variable view into a vector without copying values, and then pass
    # the flat vector to Trixi.jl's Trixi.wrap_array, which converts it into the
    # multidimensional layout expected by the spatial operators.
    return unsafe_wrap_block(evolved_variable, mesh, equations, dg, cache.cache_base)
end

# Dispatch here when semi.cache_parabolic is nothing, so the implicit wrapper has no
# HydroTrixi CacheParabolic1D (e.g. for hyperbolic or multidimensional parabolic systems).
# Trixi.default_rhs(semi_base) still uses semi_base.cache and, for parabolic systems,
# semi_base.cache_parabolic.
@inline function rhs_spatial!(du_ode, u_ode, semi_base, ::Nothing, t)
    # Preserve the blocks while the native RHS wraps and uses these flat vector aliases.
    GC.@preserve du_ode u_ode begin
        du = unsafe_vector_alias(du_ode)
        u = unsafe_vector_alias(u_ode)
        return Trixi.default_rhs(semi_base)(du, u, semi_base, t)
    end
end

function rhs_spatial!(du_ode, u_ode, semi_base, cache_parabolic::CacheParabolic1D, t)
    (; mesh, equations, boundary_conditions, source_terms, solver, solver_parabolic,
       cache) = semi_base
    GC.@preserve du_ode u_ode begin
        u = unsafe_wrap_block(u_ode, mesh, equations, solver, cache)
        du = unsafe_wrap_block(du_ode, mesh, equations, solver, cache)
        backend = Trixi.trixi_backend(u)
        time_start = time_ns()
        Trixi.@trixi_timeit Trixi.timer() "parabolic rhs!" begin
            Trixi.rhs_parabolic!(backend, du, u, t, mesh, equations, boundary_conditions,
                                  source_terms, solver, solver_parabolic, cache,
                                  cache_parabolic)
        end
        put!(semi_base.performance_counter, time_ns() - time_start)
    end
    return nothing
end

# Default operator hooks correspond to the standard semi-discrete form ∂_t u = ℛ(u, t)
@inline function nvariables_total(::AbstractTemporalOperator,
                                  semi_base)
    return Trixi.nvariables(semi_base)
end

# Fallback for the default implicit right-hand side operator (just uses the spatial RHS)
@inline function rhs_implicit!(du_ode, u_ode, ::AbstractTemporalOperator,
                               semi_base, cache_parabolic, t)
    return rhs_spatial!(du_ode, u_ode, semi_base, cache_parabolic, t)
end

# The DAE mass matrix multiplies the time derivative in M ẏ = F(y, t).
# The spatial discretization's mass matrix is handled by semi_base in rhs_spatial!.
@inline function dae_mass_matrix(u_ode, ::AbstractTemporalOperator,
                                 semi_base)
    return Diagonal(ones(eltype(u_ode), length(u_ode)))
end

function rhs_implicit!(du_ode, u_ode, operator_temporal::TemporalOperatorCapacity,
                       semi_base, cache_parabolic, t)
    rhs_spatial!(du_ode, u_ode, semi_base, cache_parabolic, t)
    (; equations) = semi_base
    capacity_function = operator_temporal.capacity_function

    # Divide by the nodal capacity after assembling the spatial residual
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
                       semi_base, cache_parabolic, t)
    evolved_variable = evolved_variable_block(u_ode, operator_temporal)
    state_variable = state_variable_block(u_ode, operator_temporal)
    evolved_variable_rhs = evolved_variable_block(du_ode, operator_temporal)
    state_variable_rhs = state_variable_block(du_ode, operator_temporal)

    # First block of du_ode: R(u_state, t)
    rhs_spatial!(evolved_variable_rhs, state_variable, semi_base, cache_parabolic, t)
    (; equations) = semi_base

    # Second block of du_ode: u_evolved - state_to_evolved(u_state)
    state_to_evolved = operator_temporal.state_to_evolved
    @inbounds for i in eachindex(state_variable_rhs, evolved_variable, state_variable)
        state_variable_rhs[i] = evolved_variable[i] -
                                state_to_evolved(state_variable[i], equations)
    end

    return nothing
end

# The DAE mass matrix gives evolved variables unit entries and algebraic state variables
# zero entries; it contains no spatial quadrature weights or element geometry factors.
function dae_mass_matrix(u_ode, ::TemporalOperatorConstitutive,
                         semi_base)
    half = length(u_ode) ÷ 2
    diagonal_entries = zeros(eltype(u_ode), length(u_ode))
    @inbounds diagonal_entries[1:half] .= one(eltype(u_ode))
    return Diagonal(diagonal_entries)
end

# SFOM variables are initialized to zero.
function passive_initial_values(::PassiveVariablesBoundaryFlux1D, semi_base, t, RealT)
    return zeros(RealT, 2)
end

function rhs_passive!(du_passive, u_physical, du_physical, ::NoPassiveVariables,
                      semi, t)
    return nothing
end

# Compute the RHS for passive variables associated with boundary fluxes (SFOM)
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

function dae_mass_matrix(u_ode, semi::SemidiscretizationImplicit)
    u_physical = physical_variable_view(u_ode, semi.passive_variables)
    physical_dae_mass_matrix = dae_mass_matrix(u_physical, semi.operator_temporal,
                                              semi.semi_base)
    n_passive = passive_variable_count(semi.passive_variables)
    if n_passive == 0
        return physical_dae_mass_matrix
    end

    # Passive scalar variables are differential variables, so we append an identity block
    passive_diagonal = ones(eltype(u_ode), n_passive)
    return Diagonal(vcat(physical_dae_mass_matrix.diag, passive_diagonal))
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
function compute_physical_coefficients(t, semi_base,
                                       ::Union{TemporalOperatorStandard,
                                               TemporalOperatorCapacity})
    return Trixi.compute_coefficients(t, semi_base)
end

function compute_physical_coefficients(t, semi_base,
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
function compute_physical_coefficients!(u_physical, t, semi_base,
                                        ::Union{TemporalOperatorStandard,
                                                TemporalOperatorCapacity})
    GC.@preserve u_physical begin
        return Trixi.compute_coefficients!(unsafe_vector_alias(u_physical), t, semi_base)
    end
end

# Compute the physical coefficients, meaning ones that are not appended passive variables
function compute_physical_coefficients!(u_physical, t,
                                        semi_base,
                                        operator_temporal::TemporalOperatorConstitutive)
    evolved_variable = evolved_variable_block(u_physical, operator_temporal)
    state_variable = state_variable_block(u_physical, operator_temporal)
    GC.@preserve u_physical begin
        Trixi.compute_coefficients!(unsafe_vector_alias(state_variable), t, semi_base)
    end
    equations = semi_base.equations
    state_to_evolved = operator_temporal.state_to_evolved

    @inbounds for i in eachindex(evolved_variable, state_variable)
        evolved_variable[i] = state_to_evolved(state_variable[i], equations)
    end

    return nothing
end

# Call compute_physical_coefficients then vcat with the passive variables
function Trixi.compute_coefficients(t, semi::SemidiscretizationImplicit)
    coefficients_physical = compute_physical_coefficients(t, semi.semi_base,
                                                          semi.operator_temporal)
    if passive_variable_count(semi.passive_variables) == 0
        return coefficients_physical
    end

    # Passive initial values are appended after the physical DAE state
    coefficients_passive = passive_initial_values(semi.passive_variables, semi.semi_base, t,
                                                  eltype(coefficients_physical))
    coefficients_ode = vcat(coefficients_physical, coefficients_passive)
    return coefficients_ode
end

function Trixi.compute_coefficients!(u_ode, t, semi::SemidiscretizationImplicit)
    u_physical = physical_variable_view(u_ode, semi.passive_variables)
    compute_physical_coefficients!(u_physical, t, semi.semi_base, semi.operator_temporal)

    if passive_variable_count(semi.passive_variables) > 0
        # Passive initial values are only written when a tail block exists
        passive_values = passive_initial_values(semi.passive_variables, semi.semi_base,
                                                t, eltype(u_ode))
        passive_variable_view(u_ode, semi) .= passive_values
    end
    return nothing
end

function rhs_implicit!(du_ode, u_ode, semi::SemidiscretizationImplicit, t)
    u_physical = physical_variable_view(u_ode, semi.passive_variables)
    du_physical = physical_variable_view(du_ode, semi.passive_variables)
    du_passive = passive_variable_view(du_ode, semi)

    # The physical residual is independent of the passive diagnostic variables
    rhs_implicit!(du_physical, u_physical, semi.operator_temporal, semi.semi_base,
                   semi.cache_parabolic, t)

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
    `SparseJacobian()` supports only serial, nonperiodic, scalar, one-dimensional `TreeMesh`
    problems using a Lobatto-Legendre `DGSEM` and `ParabolicFormulationLocalDG` with any
    penalty parameter. The supported temporal operators are
    [`TemporalOperatorStandard`](@ref), [`TemporalOperatorCapacity`](@ref), and
    [`TemporalOperatorConstitutive`](@ref); the supported passive variable types are
    [`NoPassiveVariables`](@ref) and [`PassiveVariablesBoundaryFlux1D`](@ref). MPI execution
    and periodic meshes are rejected explicitly; other configurations are unsupported and
    fail through ordinary Julia dispatch.
"""
function Trixi.semidiscretize(semi::SemidiscretizationImplicit, tspan; reset_threads = true,
                              jacobian = SparseJacobian())
    if reset_threads
        Trixi.Polyester.reset_threads!()
    end

    u0_ode = Trixi.compute_coefficients(first(tspan), semi)

    dae_mass_matrix_implicit = dae_mass_matrix(u0_ode, semi)
    jacobian_options_implicit = jacobian_options(jacobian, u0_ode, semi)
    rhs_implicit_cached = RHSImplicitCache()
    ode_function_type = SciMLBase.ODEFunction{true, SciMLBase.FullSpecialize}
    ode_function = ode_function_type(rhs_implicit_cached;
                                     mass_matrix = dae_mass_matrix_implicit,
                                     jacobian_options_implicit...)
    return SciMLBase.ODEProblem{true, SciMLBase.FullSpecialize}(ode_function, u0_ode, tspan,
                                                                semi)
end
