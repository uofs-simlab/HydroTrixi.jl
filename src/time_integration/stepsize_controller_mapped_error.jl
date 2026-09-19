# Wrap an adaptive controller to estimate error in the selected variables.
struct StepsizeControllerMappedError{Controller, ErrorControlVariables} <:
       OrdinaryDiffEqCore.AbstractController
    controller             ::Controller
    error_control_variables::ErrorControlVariables
end

function StepsizeControllerMappedError(controller, algorithm, ode,
                                       error_control_variables; reltol)
    reltol isa Number ||
        throw(ArgumentError("Mapped error control requires scalar reltol for the " *
                            "Rosenbrock linear solves."))
    algorithm isa Union{OrdinaryDiffEqRosenbrock.Rodas4,
                        OrdinaryDiffEqRosenbrock.Rodas42,
                        OrdinaryDiffEqRosenbrock.Rodas4P,
                        OrdinaryDiffEqRosenbrock.Rodas4P2,
                        OrdinaryDiffEqRosenbrock.Rodas5,
                        OrdinaryDiffEqRosenbrock.Rodas5P,
                        OrdinaryDiffEqRosenbrock.Rodas5Pe,
                        OrdinaryDiffEqRosenbrock.Rodas6P} ||
        throw(ArgumentError("error_control_variables supports Rodas4, Rodas42, " *
                            "Rodas4P, Rodas4P2, Rodas5, Rodas5P, Rodas5Pe, and Rodas6P."))
    SciMLBase.isinplace(ode) ||
        throw(ArgumentError("error_control_variables requires an in-place ODE problem."))
    if controller === nothing
        controller_type = typeof(float(first(ode.tspan)))
        controller = OrdinaryDiffEqCore.PIController(controller_type, algorithm)
    end
    return StepsizeControllerMappedError(controller, error_control_variables)
end

struct StepsizeControllerMappedErrorCache{Controller, ControllerCache,
                                         ErrorControlVariables} <:
       OrdinaryDiffEqCore.AbstractControllerCache
    # Retain the resolved controller options used by OrdinaryDiffEq.
    controller             ::Controller
    controller_cache       ::ControllerCache
    error_control_variables::ErrorControlVariables
end

function OrdinaryDiffEqCore.setup_controller_cache(algorithm, algorithm_cache,
                                                   controller::StepsizeControllerMappedError,
                                                   ::Type{ErrorEstimate},
                                                   discontinuity_problems) where {ErrorEstimate}
    controller_cache = OrdinaryDiffEqCore.setup_controller_cache(algorithm,
                                                                  algorithm_cache,
                                                                  controller.controller,
                                                                  ErrorEstimate,
                                                                  discontinuity_problems)
    return StepsizeControllerMappedErrorCache(controller_cache.controller, controller_cache,
                                              controller.error_control_variables)
end

function OrdinaryDiffEqCore.get_EEst(cache::StepsizeControllerMappedErrorCache)
    return OrdinaryDiffEqCore.get_EEst(cache.controller_cache)
end

function OrdinaryDiffEqCore.set_EEst!(cache::StepsizeControllerMappedErrorCache, error)
    return OrdinaryDiffEqCore.set_EEst!(cache.controller_cache, error)
end

function OrdinaryDiffEqCore.stepsize_controller!(integrator,
                                                 cache::StepsizeControllerMappedErrorCache,
                                                 algorithm)
    integrator.opts.step_limiter! === OrdinaryDiffEqCore.trivial_limiter! ||
        throw(ArgumentError("error_control_variables requires the default step limiter."))

    # The supported Rodas methods retain their stage increments until the controller
    # runs. Since these increments already include the time-step scaling,
    #     u - u_embedded = sum(btilde[i] * ks[i]).
    error_weights = integrator.cache.tab.btilde
    u_embedded = zero(integrator.u)
    for i in eachindex(error_weights)
        if !iszero(error_weights[i])
            @. u_embedded += error_weights[i] * integrator.cache.ks[i]
        end
    end
    @. u_embedded = integrator.u - u_embedded

    semi = integrator.p
    _, equations, _, _ = Trixi.mesh_equations_solver_cache(semi)
    u = state_variable_view(integrator.u, semi)
    u_reference = state_variable_view(u_embedded, semi)
    u_previous = state_variable_view(integrator.uprev, semi)
    # Array tolerances follow the full ODE storage layout, as required by SciML.
    abstol = integrator.opts.abstol
    if !(abstol isa Number)
        abstol = state_variable_view(abstol, semi)
    end

    mapping = cache.error_control_variables
    reltol = integrator.opts.reltol
    residual = broadcast(u, u_reference, u_previous,
                         abstol, reltol) do value, reference, previous, atol, rtol
        mapped_value = mapping(value, equations)
        mapped_reference = mapping(reference, equations)
        mapped_previous = mapping(previous, equations)
        return (mapped_value .- mapped_reference) ./
               (atol .+ rtol .* max.(abs.(mapped_previous), abs.(mapped_value)))
    end
    error = Trixi.ode_norm(residual, integrator.t + integrator.dt)
    OrdinaryDiffEqCore.set_EEst!(integrator, error)
    return OrdinaryDiffEqCore.stepsize_controller!(integrator, cache.controller_cache,
                                                   algorithm)
end

function OrdinaryDiffEqCore.accept_step_controller(integrator,
                                                   cache::StepsizeControllerMappedErrorCache,
                                                   algorithm)
    return OrdinaryDiffEqCore.accept_step_controller(integrator, cache.controller_cache,
                                                     algorithm)
end

function OrdinaryDiffEqCore.step_accept_controller!(integrator,
                                                    cache::StepsizeControllerMappedErrorCache,
                                                    algorithm, q)
    return OrdinaryDiffEqCore.step_accept_controller!(integrator, cache.controller_cache,
                                                      algorithm, q)
end

function OrdinaryDiffEqCore.step_reject_controller!(integrator,
                                                    cache::StepsizeControllerMappedErrorCache,
                                                    algorithm)
    return OrdinaryDiffEqCore.step_reject_controller!(integrator, cache.controller_cache,
                                                      algorithm)
end

function OrdinaryDiffEqCore.reinit_controller!(integrator::SciMLBase.DEIntegrator,
                                               cache::StepsizeControllerMappedErrorCache)
    return OrdinaryDiffEqCore.reinit_controller!(integrator, cache.controller_cache)
end
