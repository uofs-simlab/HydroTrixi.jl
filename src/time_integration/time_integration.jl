@doc raw"""
    default_algorithm(ode::SciMLBase.ODEProblem; kwargs...)

Return a recommended OrdinaryDiffEq.jl time integration algorithm for `ode`, suitable for
passing to `SciMLBase.solve`.

The recommended algorithm for `SemidiscretizationImplicit` is `Rodas5P`, an eight-stage,
fifth-order Rosenbrock-Wanner method, and it recomputes the Jacobian after at most one time
step. When `ode` has the sparse Jacobian prototype supplied by [`semidiscretize`](@ref) with
[`SparseJacobian`](@ref), the algorithm uses sparse forward-mode automatic differentiation
and a KLU linear solver. Otherwise, it uses dense forward-mode automatic differentiation
and a dense LU linear solver.

Keyword arguments override these defaults or are forwarded to the `Rodas5P` constructor.

# References
- Steinebach, G. (2023). Construction of Rosenbrock-Wanner method Rodas5P and numerical
  benchmarks within the Julia Differential Equations package. *BIT Numerical
  Mathematics*, 63, Article 27.
  [DOI: 10.1007/s10543-023-00967-x](https://doi.org/10.1007/s10543-023-00967-x)
"""
function default_algorithm(ode::SciMLBase.ODEProblem{U, T, I, P};
                           chunk_size = Val{0}(),
                           autodiff = OrdinaryDiffEqRosenbrock.AutoForwardDiff(),
                           standardtag = Val{true}(),
                           diff_type = Val{:forward}(),
                           linsolve = ode.f.jac_prototype isa SparseMatrixCSC ?
                                      LinearSolve.KLUFactorization() :
                                      LinearSolve.LUFactorization(),
                           precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
                           step_limiter! = OrdinaryDiffEqCore.trivial_limiter!,
                           stage_limiter! = OrdinaryDiffEqCore.trivial_limiter!,
                           concrete_jac = nothing,
                           max_jac_age = 1,
                           jac_reuse_gamma_tol = 0.03,
                           kwargs...) where {U, T, I, P <: SemidiscretizationImplicit}
    return OrdinaryDiffEqRosenbrock.Rodas5P(; chunk_size, autodiff, standardtag,
                                            concrete_jac, diff_type, linsolve, precs,
                                            step_limiter!, stage_limiter!, max_jac_age,
                                            jac_reuse_gamma_tol, kwargs...)
end

@doc raw"""
    state_variable_norm(semi::SemidiscretizationImplicit)

Return a norm that restricts adaptive error control to the state-variable degrees of
freedom in `semi`. The returned callable can be passed as the `internalnorm` keyword to
the SciML `solve` function as follows:
```julia
error_norm = state_variable_norm(semi)
sol = solve(ode, default_algorithm(ode); internalnorm = error_norm, kwargs...)
```

For [`TemporalOperatorStandard`](@ref) and [`TemporalOperatorCapacity`](@ref), the norm
uses the complete physical state and excludes appended passive diagnostic variables. For
[`TemporalOperatorConstitutive`](@ref), it uses only the state-variable block and excludes
the evolved-variable block and passive diagnostic variables. Thus, for the Richards
equation, it restricts error control to pressure head for both [`PressureHeadForm`](@ref)
and [`MixedForm`](@ref).
"""
function state_variable_norm(semi::SemidiscretizationImplicit)
    return function (u, t)
        u isa Number && return Trixi.ode_norm(u, t)
        state_variables = state_variable_view(u, semi)
        return Trixi.ode_norm(state_variables, t)
    end
end
