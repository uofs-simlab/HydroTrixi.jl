@inline function Trixi.analyze(::typeof(water_content), du::AbstractVector,
                               u_ode::AbstractVector, t,
                               semi::SemidiscretizationImplicit)
    return first(evolved_variables_integral(u_ode, semi))
end

@inline function Trixi.analyze(::typeof(water_content), du, u, t,
                               semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return first(Trixi.integrate(Trixi.cons2cons, u, mesh, equations, solver, cache;
                                 normalize = false))
end

Trixi.pretty_form_ascii(::typeof(water_content)) = "water_content"
Trixi.pretty_form_utf(::typeof(water_content)) = "∫θ"

function water_content_timederivative end

@inline function Trixi.analyze(::typeof(water_content_timederivative),
                               du_ode::AbstractVector,
                               u_ode::AbstractVector, t,
                               semi::SemidiscretizationImplicit)
    return first(Trixi.integrate(Trixi.cons2cons, du_ode, semi; normalize = false))
end

@inline function Trixi.analyze(::typeof(water_content_timederivative), du, u, t,
                               semi::SemidiscretizationImplicit)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi.semi_base)
    return first(Trixi.integrate(Trixi.cons2cons, du, mesh, equations, solver, cache;
                                 normalize = false))
end

Trixi.pretty_form_ascii(::typeof(water_content_timederivative)) = "water_content_t"
Trixi.pretty_form_utf(::typeof(water_content_timederivative)) = "d/dt ∫θ"
