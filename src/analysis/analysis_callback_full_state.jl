@doc raw"""
    AnalysisCallbackFullState(ode; full_state_analysis_integrals, kwargs...)

Create a Trixi analysis callback that additionally supports scalar analysis quantities
requiring the complete ODE state. The remaining keyword arguments are forwarded to
`Trixi.AnalysisCallback`.

Before evaluating its ordinary analysis integrals, `Trixi.AnalysisCallback` calls
`Trixi.wrap_array` once for the complete ODE vector. HydroTrixi specializes that operation
for `SemidiscretizationImplicit`: standard and capacity formulations expose the physical
state after removing appended passive variables, while constitutive formulations expose
only the evolved-variable block after also removing the state-variable block. Information
discarded by this shared wrapped view cannot be recovered by an individual ordinary
analysis integral.

This callback therefore evaluates each quantity in `full_state_analysis_integrals` as

```julia
Trixi.analyze(quantity, du_ode, u_ode, t, semi)
```

where `du_ode` and `u_ode` are the complete solver vectors, including any algebraic or
passive variables appended by the semidiscretization. Ordinary analysis errors and
integrals continue to use `Trixi.AnalysisCallback` unchanged.

The callback owns only adapters for routing the complete state. The numerical arrays
remain owned by the integrator and are refreshed in the adapters immediately before each
analysis invocation.
"""
struct AnalysisCallbackFullState{Callback, FullStateAnalysisIntegrals}
    callback::Callback
    full_state_analysis_integrals::FullStateAnalysisIntegrals
end

# Adapt a full-state quantity to Trixi's standard analysis-integral calling convention.
mutable struct FullStateAnalysisIntegral{Quantity, Du, U}
    quantity::Quantity
    du_ode::Du
    u_ode::U
end

function FullStateAnalysisIntegral(quantity, u_ode)
    return FullStateAnalysisIntegral(quantity, similar(u_ode), u_ode)
end

function Base.show(io::IO, integral::FullStateAnalysisIntegral)
    return show(io, integral.quantity)
end

function Trixi.pretty_form_ascii(integral::FullStateAnalysisIntegral)
    return Trixi.pretty_form_ascii(integral.quantity)
end

function Trixi.pretty_form_utf(integral::FullStateAnalysisIntegral)
    return Trixi.pretty_form_utf(integral.quantity)
end

@inline function Trixi.analyze(integral::FullStateAnalysisIntegral, du, u, t,
                               semi::Trixi.AbstractSemidiscretization)
    return Trixi.analyze(integral.quantity, integral.du_ode, integral.u_ode, t, semi)
end

@inline function update_full_state!(integral::FullStateAnalysisIntegral, du_ode, u_ode)
    integral.du_ode = du_ode
    integral.u_ode = u_ode
    return nothing
end

@inline update_full_state!(::Tuple{}, du_ode, u_ode) = nothing

@inline function update_full_state!(integrals::Tuple, du_ode, u_ode)
    update_full_state!(first(integrals), du_ode, u_ode)
    update_full_state!(Base.tail(integrals), du_ode, u_ode)
    return nothing
end

function AnalysisCallbackFullState(ode::SciMLBase.ODEProblem;
                                   full_state_analysis_integrals,
                                   extra_analysis_integrals = (), kwargs...)
    u_ode = ode.u0
    full_state_integrals = map(quantity -> FullStateAnalysisIntegral(quantity, u_ode),
                               Tuple(full_state_analysis_integrals))
    analysis_integrals = (Tuple(extra_analysis_integrals)..., full_state_integrals...)

    # Create the Trixi analysis callback with the standard analysis integrals
    callback = Trixi.AnalysisCallback(ode.p;
                                      extra_analysis_integrals = analysis_integrals,
                                      kwargs...)

    # Create the extension of the analysis callback to use full-state analysis integrals
    analysis_callback = AnalysisCallbackFullState(callback, full_state_integrals)

    return SciMLBase.DiscreteCallback(callback.condition, analysis_callback;
                                      initialize = initialize_analysis_callback_full_state!,
                                      finalize = finalize_analysis_callback_full_state!,
                                      save_positions = callback.save_positions,
                                      initializealg = callback.initializealg,
                                      saved_clock_partitions = callback.saved_clock_partitions,
                                      initialize_save_discretes =
                                      callback.initialize_save_discretes)
end

# Refresh complete-state references before delegating to Trixi's callback implementation.
function (analysis_callback::AnalysisCallbackFullState)(integrator)
    u_ode = integrator.u
    du_ode = first(SciMLBase.get_tmp_cache(integrator))
    update_full_state!(analysis_callback.full_state_analysis_integrals, du_ode, u_ode)
    return analysis_callback.callback.affect!(integrator)
end

# Preserve Trixi's postprocessing interface used for error checks and EOC analysis.
function (cb::SciMLBase.DiscreteCallback{Condition, Affect!})(sol) where {
        Condition, Affect! <: AnalysisCallbackFullState}
    return cb.affect!.callback(sol)
end

function initialize_analysis_callback_full_state!(cb, u_ode, t, integrator)
    analysis_callback = cb.affect!
    du_ode = first(SciMLBase.get_tmp_cache(integrator))
    update_full_state!(analysis_callback.full_state_analysis_integrals, du_ode, u_ode)
    callback = analysis_callback.callback
    return callback.initialize(callback, u_ode, t, integrator)
end

function finalize_analysis_callback_full_state!(cb, u_ode, t, integrator)
    callback = cb.affect!.callback
    return callback.finalize(callback, u_ode, t, integrator)
end

function Base.show(io::IO, mime::MIME"text/plain",
                   cb::SciMLBase.DiscreteCallback{<:Any,
                                                  <:AnalysisCallbackFullState})
    return show(io, mime, cb.affect!.callback)
end
