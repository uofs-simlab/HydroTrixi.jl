# The nodal water-content vector θ_k is stored directly in the mixed formulation and
# evaluated as ϑ(ψ_k) in the pressure-head formulation.
@inline function water_content_integral(u_ode::AbstractVector,
                                        semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)
    return water_content_integral(u, mesh, equations, solver, cache.cache_base,
                                  semi.operator_temporal)
end

# Destructure the semidiscretization so the integral can dispatch on the equations and
# temporal operator.
@inline function water_content_integral(u, semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    return water_content_integral(u, mesh, equations, solver, cache.cache_base,
                                  semi.operator_temporal)
end

# In mixed form, Trixi.wrap_array extracts the evolved water-content vector θ_k.
@inline function water_content_integral(u, mesh, equations::RichardsEquation1D, solver,
                                        cache,
                                        ::TemporalOperatorConstitutive)
    return first(Trixi.integrate(Trixi.cons2cons, u, mesh, equations, solver, cache;
                                 normalize = false))
end

# In pressure-head form, evaluate ϑ(ψ_k) at the pressure-head nodes.
@inline function water_content_integral(u, mesh, equations::RichardsEquation1D, solver,
                                        cache,
                                        ::TemporalOperatorCapacity)
    return Trixi.integrate(water_content, u, mesh, equations, solver, cache;
                           normalize = false)
end

@inline function Trixi.analyze(::typeof(water_content), du, u, t,
                               semi::SemidiscretizationImplicit)
    return water_content_integral(u, semi)
end

Trixi.pretty_form_ascii(::typeof(water_content)) = "water_content"
Trixi.pretty_form_utf(::typeof(water_content)) = "∫θ"

@doc raw"""
    water_content_timederivative

Analysis integral for the time derivative of the total water content,
```math
\frac{\mathrm{d}}{\mathrm{d}t}
\sum_{k=1}^{K}J_k\boldsymbol{1}^{\mathrm{T}}\boldsymbol{W}
\boldsymbol{\theta}_k(t).
```

For the mixed formulation, ``\boldsymbol{\theta}_k`` is evolved directly. For the
pressure-head formulation,
``\boldsymbol{\theta}_k=\boldsymbol{\vartheta}(\boldsymbol{\psi}_k)`` and this quantity
integrates ``\boldsymbol{C}_k\dot{\boldsymbol{\psi}}_k``, where
``\boldsymbol{C}_k`` contains the nodal values of [`water_capacity`](@ref). The two
expressions agree when the constitutive constraint
``\boldsymbol{\theta}_k=\boldsymbol{\vartheta}(\boldsymbol{\psi}_k)`` is satisfied.
"""
function water_content_timederivative end

# The manuscript denotes the nodal rate by θ̇_k in mixed form and C_k ψ̇_k in
# pressure-head form.
@inline function water_content_timederivative_integral(du_ode::AbstractVector,
                                                       u_ode::AbstractVector,
                                                       semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)
    du = Trixi.wrap_array(du_ode, mesh, equations, solver, cache)
    return water_content_timederivative_integral(du, u, mesh, equations, solver,
                                                 cache.cache_base,
                                                 semi.operator_temporal)
end

@inline function water_content_timederivative_integral(du, u,
                                                       semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    return water_content_timederivative_integral(du, u, mesh, equations, solver,
                                                 cache.cache_base,
                                                 semi.operator_temporal)
end

@inline function water_content_timederivative_integral(du, u, mesh,
                                                       equations::RichardsEquation1D,
                                                       solver, cache,
                                                       ::TemporalOperatorConstitutive)
    return first(Trixi.integrate(Trixi.cons2cons, du, mesh, equations, solver, cache;
                                 normalize = false))
end

@inline function water_content_timederivative_integral(du, u, mesh,
                                                       equations::RichardsEquation1D,
                                                       solver, cache,
                                                       ::TemporalOperatorCapacity)
    return Trixi.integrate_via_indices(u, mesh, equations, solver, cache, du;
                                       normalize = false) do u, i, element, equations,
                                                             solver, du
        u_node = Trixi.get_node_vars(u, equations, solver, i, element)
        du_node = Trixi.get_node_vars(du, equations, solver, i, element)
        return water_capacity(u_node, equations) * pressure_head(du_node)
    end
end

@inline function Trixi.analyze(::typeof(water_content_timederivative), du, u, t,
                               semi::SemidiscretizationImplicit)
    return water_content_timederivative_integral(du, u, semi)
end

Trixi.pretty_form_ascii(::typeof(water_content_timederivative)) = "water_content_t"
Trixi.pretty_form_utf(::typeof(water_content_timederivative)) = "d/dt ∫θ"

@doc raw"""
    mass_balance
    mass_balance(u_ode, semi::SemidiscretizationImplicit)

Return the signed cumulative water mass balance
```math
F_{\mathrm{B}}(t)-F_{\mathrm{T}}(t)
-\sum_{k=1}^{K}J_k\boldsymbol{1}^{\mathrm{T}}\boldsymbol{W}
 \boldsymbol{\theta}_k(t),
```
where ``F_{\mathrm{T}}`` and ``F_{\mathrm{B}}`` are the passive SFOM variables satisfying
``\dot{F}_{\mathrm{T}}=f_{\mathrm{T}}^\star`` and
``\dot{F}_{\mathrm{B}}=f_{\mathrm{B}}^\star``, respectively, and
``\boldsymbol{\theta}_k`` is the nodal water-content vector. The Julia named-tuple keys
`x_neg` and `x_pos` represent the top and bottom boundaries. This function returns the
negative of the fully discrete invariant
```math
\sum_{k=1}^{K}J_k\boldsymbol{1}^{\mathrm{T}}\boldsymbol{W}
\boldsymbol{\theta}_k(t)+F_{\mathrm{T}}(t)-F_{\mathrm{B}}(t)
```
used in the accompanying manuscript. Differencing this quantity between time levels
produces the signed mass bias ``\epsilon_{\mathrm{B}}^n`` documented by
[`mass_bias`](@ref).

This solver flux output method diagnostic requires
[`PassiveVariablesBoundaryFlux1D`](@ref), which advances the integrated numerical
boundary fluxes with the same time-integration stages, time steps, and step-acceptance
decisions as the physical state.
"""
function mass_balance end

function mass_balance(u_ode::AbstractVector, semi::SemidiscretizationImplicit)
    boundary_fluxes = boundary_flux_integrals(u_ode, semi)
    storage = water_content_integral(u_ode, semi)
    return -boundary_fluxes.x_neg + boundary_fluxes.x_pos - storage
end

@inline function Trixi.analyze(::typeof(mass_balance), du_ode::AbstractVector,
                               u_ode::AbstractVector, t,
                               semi::SemidiscretizationImplicit)
    return mass_balance(u_ode, semi)
end

Trixi.pretty_form_ascii(::typeof(mass_balance)) = "mass_balance"
Trixi.pretty_form_utf(::typeof(mass_balance)) = "mass balance"

@doc raw"""
    mass_bias(u_ode, semi::SemidiscretizationImplicit, initial_water_content)

Return the water mass bias relative to an initial total water content:
```math
\epsilon_{\mathrm{B}}^n =
\left(F_{\mathrm{B}}^n-F_{\mathrm{B}}^0\right)
-\left(F_{\mathrm{T}}^n-F_{\mathrm{T}}^0\right)
-\left(
\sum_{k=1}^{K^n}J_k^n\boldsymbol{1}^{\mathrm{T}}\boldsymbol{W}
\boldsymbol{\theta}_k^n
-
\sum_{k=1}^{K^0}J_k^0\boldsymbol{1}^{\mathrm{T}}\boldsymbol{W}
\boldsymbol{\theta}_k^0
\right).
```
This method assumes the initialization
``F_{\mathrm{T}}^0=F_{\mathrm{B}}^0=0``. For complete
saved time histories, prefer [`mass_bias_history`](@ref), which differences
[`mass_balance`](@ref) and therefore also supports nonzero initial SFOM variables.
"""
function mass_bias end

function mass_bias(u_ode::AbstractVector, semi::SemidiscretizationImplicit,
                   initial_water_content)
    return mass_balance(u_ode, semi) + initial_water_content
end

@doc raw"""
    mass_bias_history(sol; initial_water_content = nothing)
    mass_bias_history(analysis_path::AbstractString;
                      time_column = "time", mass_balance_column = "mass_balance")

Return the saved times and corresponding water mass biases for `sol`.

The solution must use a [`SemidiscretizationImplicit`](@ref) with
[`PassiveVariablesBoundaryFlux1D`](@ref). By default, the first saved state must be at the
initial time ``t^0``. When it is not available, pass `initial_water_content` explicitly;
this assumes ``F_{\mathrm{T}}^0=F_{\mathrm{B}}^0=0``. Solution
postprocessing requires a fixed state layout; for AMR solutions, use the analysis-file
method instead.

When `analysis_path` is provided, read the time and signed cumulative water mass balance
columns from an analysis file written by
`AnalysisCallbackFullState(save_analysis = true)` and subtract the first value.
"""
function mass_bias_history(sol; initial_water_content = nothing)
    semi = sol.prob.p
    state_length = length(last(sol.u))
    if any(u_ode -> length(u_ode) != state_length, sol.u)
        throw(ArgumentError("Saved states have different sizes. Use the analysis-file " *
                            "method for AMR solutions."))
    end

    if !isnothing(initial_water_content)
        biases = [mass_bias(u_ode, semi, initial_water_content) for u_ode in sol.u]
        return collect(sol.t), biases
    end

    if first(sol.t) != first(sol.prob.tspan)
        throw(ArgumentError("The solution does not contain the initial state. Pass " *
                            "`initial_water_content` explicitly."))
    end

    balances = [mass_balance(u_ode, semi) for u_ode in sol.u]
    biases = balances .- first(balances)

    return collect(sol.t), biases
end

function mass_bias_history(analysis_path::AbstractString; time_column = "time",
                           mass_balance_column = "mass_balance")
    times = Float64[]
    balances = Float64[]

    open(analysis_path, "r") do io
        # Read the Trixi.jl analysis header to locate scalar output columns
        header = nothing
        for line in eachline(io)
            stripped_line = strip(line)
            if isempty(stripped_line)
                continue
            end
            header = stripped_line
            break
        end

        header_columns = split(strip(header[2:end]))
        column_indices = Dict(column => index for (index, column) in pairs(header_columns))
        time_index = column_indices[time_column]
        mass_balance_index = column_indices[mass_balance_column]

        # Parse the scalar time history from the selected columns
        for line in eachline(io)
            stripped_line = strip(line)
            if isempty(stripped_line)
                continue
            end
            if startswith(stripped_line, "#")
                continue
            end

            values = split(stripped_line)
            push!(times, parse(Float64, values[time_index]))
            push!(balances, parse(Float64, values[mass_balance_index]))
        end
    end

    if isempty(times)
        throw(ArgumentError("`analysis_path` does not contain mass-balance samples."))
    end

    biases = balances .- first(balances)
    return times, biases
end
