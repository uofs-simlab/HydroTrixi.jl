# Wrap an adaptive controller to estimate error in the selected variables.
struct StepsizeControllerMappedError{Controller, ErrorControlBlock, ErrorControlMapping} <:
       OrdinaryDiffEqCore.AbstractController
    controller           ::Controller
    error_control_block  ::ErrorControlBlock
    error_control_mapping::ErrorControlMapping
end

function StepsizeControllerMappedError(controller, algorithm, ode,
                                       error_control_block, error_control_mapping; reltol)
    if !(reltol isa Number)
        throw(ArgumentError("Mapped error control requires scalar reltol for the " *
                            "Rosenbrock linear solves."))
    end
    if !(algorithm isa Union{OrdinaryDiffEqRosenbrock.Rodas4,
                             OrdinaryDiffEqRosenbrock.Rodas42,
                             OrdinaryDiffEqRosenbrock.Rodas4P,
                             OrdinaryDiffEqRosenbrock.Rodas4P2,
                             OrdinaryDiffEqRosenbrock.Rodas5,
                             OrdinaryDiffEqRosenbrock.Rodas5P,
                             OrdinaryDiffEqRosenbrock.Rodas5Pe,
                             OrdinaryDiffEqRosenbrock.Rodas6P})
        throw(ArgumentError("error_control_mapping supports Rodas4, Rodas42, " *
                            "Rodas4P, Rodas4P2, Rodas5, Rodas5P, Rodas5Pe, and Rodas6P."))
    end
    if !SciMLBase.isinplace(ode)
        throw(ArgumentError("error_control_mapping requires an in-place ODE problem."))
    end
    if controller === nothing
        controller_type = typeof(float(first(ode.tspan)))
        controller = OrdinaryDiffEqCore.PIController(controller_type, algorithm)
    end
    return StepsizeControllerMappedError(controller, error_control_block,
                                         error_control_mapping)
end

struct StepsizeControllerMappedErrorCache{Controller, ControllerCache, ErrorControlBlock,
                                          ErrorControlMapping} <:
       OrdinaryDiffEqCore.AbstractControllerCache
    # Retain the resolved controller options used by OrdinaryDiffEq.
    controller           ::Controller
    controller_cache     ::ControllerCache
    error_control_block  ::ErrorControlBlock
    error_control_mapping::ErrorControlMapping
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
                                              controller.error_control_block,
                                              controller.error_control_mapping)
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
    if integrator.opts.step_limiter! !== OrdinaryDiffEqCore.trivial_limiter!
        throw(ArgumentError("error_control_mapping requires the default step limiter."))
    end

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
    block = cache.error_control_block
    u = block(integrator.u, semi)
    u_reference = block(u_embedded, semi)
    u_previous = block(integrator.uprev, semi)
    # Array tolerances follow the full ODE storage layout, as required by SciML.
    abstol = integrator.opts.abstol
    if !(abstol isa Number)
        abstol = block(abstol, semi)
    end

    mapping = cache.error_control_mapping
    reltol = integrator.opts.reltol
    residual = broadcast(u, u_reference, u_previous,
                         abstol, reltol) do value, reference, previous, atol, rtol
        mapped_value = mapping(value, equations)
        mapped_reference = mapping(reference, equations)
        mapped_previous = mapping(previous, equations)
        return (mapped_value - mapped_reference) /
               (atol + rtol * max(abs(mapped_previous), abs(mapped_value)))
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
