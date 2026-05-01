function evolved_variables_integral(u_ode, semi::SemidiscretizationImplicit)
    return Trixi.integrate(Trixi.cons2cons, u_ode, semi; normalize = false)
end

function evolved_variables_timederivative(u_ode, semi::SemidiscretizationImplicit, t)
    du_ode = similar(u_ode)
    rhs_implicit!(du_ode, u_ode, semi, t)
    return Trixi.integrate(Trixi.cons2cons, du_ode, semi; normalize = false)
end

include("analysis_richards.jl")
