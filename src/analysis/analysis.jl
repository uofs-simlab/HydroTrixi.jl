function evolved_variables_integral(u_ode, semi::SemidiscretizationImplicit)
    return Trixi.integrate(Trixi.cons2cons, u_ode, semi; normalize = false)
end

include("analysis_richards.jl")
