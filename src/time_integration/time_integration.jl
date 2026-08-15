@doc raw"""
    default_algorithm(semi::SemidiscretizationImplicit)

Return a recommended OrdinaryDiffEq.jl time integration algorithm for the implicit
semidiscretization `semi`, suitable for passing to `SciMLBase.solve`.

The recommended algorithm is `Rodas5P()`, an eight-stage, fifth-order Rosenbrock-Wanner
method with graph-coloured forward-mode automatic differentiation when
[`semidiscretize`](@ref) supplies the default [`SparseJacobian`](@ref) prototype. The Jacobian strategy passed to [`semidiscretize`](@ref) configures Jacobian storage. The
time integration algorithm configures the differentiation backend and linear solver.

# References
- Steinebach, G. (2023). Construction of Rosenbrock-Wanner method Rodas5P and numerical
  benchmarks within the Julia Differential Equations package. *BIT Numerical
  Mathematics*, 63, Article 27.
  [DOI: 10.1007/s10543-023-00967-x](https://doi.org/10.1007/s10543-023-00967-x)
"""
function default_algorithm(::SemidiscretizationImplicit)
    return OrdinaryDiffEqRosenbrock.Rodas5P()
end

@doc raw"""
    state_variable_norm(semi::SemidiscretizationImplicit)

Return a norm that restricts adaptive error control to the state-variable degrees of
freedom in `semi`. The returned callable can be passed as the `internalnorm` keyword to
the SciML `solve` function as follows:
```julia
error_norm = state_variable_norm(semi)
sol = solve(ode, default_algorithm(semi); internalnorm = error_norm, kwargs...)
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
