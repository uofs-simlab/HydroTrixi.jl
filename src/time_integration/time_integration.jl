@doc raw"""
    default_algorithm(ode::SciMLBase.ODEProblem; kwargs...)

Return a recommended OrdinaryDiffEq.jl time integration algorithm for `ode`, suitable for
passing to `SciMLBase.solve`.

The recommended algorithm for `SemidiscretizationImplicit` is `Rodas5P`, an eight-stage,
fifth-order Rosenbrock-Wanner method, and it recomputes the Jacobian after at most one time
step. When `ode` has the sparse Jacobian prototype supplied by [`semidiscretize`](@ref) with
[`SparseJacobian`](@ref), the algorithm uses sparse forward-mode automatic differentiation
with a deterministic analytical colouring and a KLU linear solver. The colouring keeps
the sparse differentiation cache type unchanged when AMR changes a mesh that is large
enough to contain the complete colour palette. Otherwise, the algorithm uses dense
forward-mode automatic differentiation with automatic chunk-size selection and a dense LU
linear solver.

Keyword arguments override these defaults or are forwarded to the `Rodas5P` constructor.
Configure differentiation through `autodiff` and preconditioning through `linsolve`.
Pass `step_limiter` and `stage_limiter` to `solve` or `solve_implicit`.

# References
- Davis, T. A., and Palamadai Natarajan, E. (2010). Algorithm 907: KLU, a direct sparse
  solver for circuit simulation problems. *ACM Transactions on Mathematical Software*,
  37(3), Article 36.
  [DOI: 10.1145/1824801.1824814](https://doi.org/10.1145/1824801.1824814)
- Revels, J., Lubin, M., and Papamarkou, T. (2016). Forward-mode automatic differentiation
  in Julia. *arXiv:1607.07892*.
  [DOI: 10.48550/arXiv.1607.07892](https://doi.org/10.48550/arXiv.1607.07892)
- Steinebach, G. (2023). Construction of Rosenbrock-Wanner method Rodas5P and numerical
  benchmarks within the Julia Differential Equations package. *BIT Numerical
  Mathematics*, 63, Article 27.
  [DOI: 10.1007/s10543-023-00967-x](https://doi.org/10.1007/s10543-023-00967-x)
"""
function default_algorithm(ode::SciMLBase.ODEProblem{U, T, I, P};
                           autodiff = OrdinaryDiffEqRosenbrock.AutoForwardDiff(),
                           linsolve = ode.f.jac_prototype isa SparseMatrixCSC ?
                                      LinearSolve.KLUFactorization() :
                                      LinearSolve.LUFactorization(),
                           concrete_jac = nothing,
                           max_jac_age = 1,
                           jac_reuse_gamma_tol = 0.03,
                           kwargs...) where {U, T, I, P <: SemidiscretizationImplicit}
    return OrdinaryDiffEqRosenbrock.Rodas5P(; autodiff, concrete_jac, linsolve, max_jac_age,
                                            jac_reuse_gamma_tol, kwargs...)
end

"""
    pressure_head_out_of_domain(u, semi, t)

Return `true` if any pressure-head degree of freedom in the candidate ODE state `u` is
nonnegative. The call signature matches SciML's `isoutofdomain` predicate and can be
passed directly to [`solve_implicit`](@ref):
```julia
sol = solve_implicit(ode; isoutofdomain = pressure_head_out_of_domain, kwargs...)
```

The domain criterion dispatches on the equations stored by
[`SemidiscretizationImplicit`](@ref). For [`RichardsEquation1D`](@ref), the pressure head
is obtained from the state-variable block, so the criterion applies to both the mixed and
pressure-head forms. Returning `true` rejects the candidate time step; it does not
constrain internal time-integration stages.
"""
@inline function pressure_head_out_of_domain(u, semi::SemidiscretizationImplicit, t)
    _, equations, _, _ = Trixi.mesh_equations_solver_cache(semi)
    return pressure_head_out_of_domain(u, semi, t, equations)
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

@doc raw"""
    evolved_variable_norm(semi::SemidiscretizationImplicit)

Return a norm that restricts adaptive error control to the evolved-variable degrees of
freedom in `semi`. The returned callable can be passed as the `internalnorm` keyword to
the SciML `solve` function as follows:
```julia
error_norm = evolved_variable_norm(semi)
sol = solve(ode, default_algorithm(ode); internalnorm = error_norm, kwargs...)
```

For [`TemporalOperatorStandard`](@ref) and [`TemporalOperatorCapacity`](@ref), the norm
uses the complete physical state and excludes appended passive diagnostic variables. For
[`TemporalOperatorConstitutive`](@ref), it uses only the evolved-variable block and
excludes the state-variable block and passive diagnostic variables. Thus, for the Richards
equation, it restricts error control to pressure head for [`PressureHeadForm`](@ref) and
water content for [`MixedForm`](@ref).
"""
function evolved_variable_norm(semi::SemidiscretizationImplicit)
    return function (u, t)
        u isa Number && return Trixi.ode_norm(u, t)
        evolved_variables = evolved_variable_view(u, semi)
        return Trixi.ode_norm(evolved_variables, t)
    end
end

"""
    default_stepsize_controller(algorithm, ode)

Return HydroTrixi.jl's default adaptive step-size controller for `algorithm` and `ode`.

For `Rodas5P`, return the PI controller used by [`solve_implicit`](@ref). For other
algorithms, return `nothing` so that OrdinaryDiffEq.jl selects the algorithm's default
controller.
"""
default_stepsize_controller(algorithm, ode) = nothing

function default_stepsize_controller(algorithm::OrdinaryDiffEqRosenbrock.Rodas5P,
                                     ode)
    controller_type = typeof(float(first(ode.tspan)))
    return OrdinaryDiffEqCore.PIController(controller_type, algorithm;
                                           # Current- and previous-error exponents
                                           beta1 = 0.14,
                                           beta2 = 0.08,
                                           # Minimum shrink and maximum growth factors
                                           qmin = 0.2,
                                           qmax = 10.0,
                                           # Allow larger growth after the first step
                                           qmax_first_step = 1.0e4,
                                           # Safety factor
                                           gamma = 0.9,
                                           # OrdinaryDiffEq's deadband, which holds the
                                           # time step fixed when the controller gives
                                           # a time-step divisor between 1.0 and 1.2
                                           qsteady_min = 1.0,
                                           qsteady_max = 1.2,
                                           # Previous-error initialization and floor
                                           qoldinit = 1.0e-4)
end

@inline function error_control_uses_evolved_variables(mapping,
                                                      operator::TemporalOperatorConstitutive)
    return mapping === operator.state_to_evolved
end

@inline error_control_uses_evolved_variables(mapping, operator) = false

@doc raw"""
    solve_implicit(ode, algorithm=default_algorithm(ode);
                   dt, adaptive=true, abstol=1.0e-11, reltol=1.0e-7,
                   error_control_variables=nothing,
                   kwargs...)

Solve `ode` with HydroTrixi.jl's implicit time integration defaults. The initial time step
`dt` is required. The absolute and relative tolerances default to `1.0e-11` and `1.0e-7`,
respectively.

The adaptive defaults use a PI controller configured for the fifth-order
[`default_algorithm`](@ref), with coefficients ``\beta_1=0.14`` and ``\beta_2=0.08``,
safety factor `0.9`, maximum growth factor `10`, maximum shrink factor `0.2`, and initial
and minimum stored previous error `1.0e-4`. The maximum growth factor is `1.0e4` for the
first step-size proposal and `10` thereafter. OrdinaryDiffEq.jl's default steady-step
deadband holds the time step fixed when the controller proposes a time-step divisor
between `1.0` and `1.2`. These parameters are specified explicitly by
[`default_stepsize_controller`](@ref).

The defaults also use [`evolved_variable_norm`](@ref), disable saving every accepted step,
and allow `typemax(Int)` iterations. Consistent initial states are constructed by
[`semidiscretize`](@ref), and AMR reconstructs state variables from the transferred evolved
variables. No additional nonlinear correction is requested; `initializealg` can be passed
through to SciML to override its initialization method.

`error_control_variables = nothing` means no additional conversion: the existing
`internalnorm` controls the error. Its default, [`evolved_variable_norm`](@ref), selects
stored water content for [`MixedForm`](@ref) and stored pressure head for
[`PressureHeadForm`](@ref). A user-supplied `internalnorm` is preserved.

With an explicit function `error_control_variables(value, equations)`, apply that
function pointwise to the physical state-variable entries of the candidate, embedded,
and previous solutions before forming and scaling the error. Passive diagnostics are
excluded. For example, on a pressure-head problem:
```julia
sol = solve_implicit(ode; dt = 1.0e-2, error_control_variables = water_content,
                     abstol = 1.0e-8, reltol = 1.0e-5)
```
The tolerances then apply to the converted quantities. Scalar tolerances are broadcast;
array `abstol` must follow the full ODE layout, and its state-variable entries are
selected. `reltol` must be scalar because it is also used by the Rosenbrock linear solves.
An explicit function always requests conversion, including for [`MixedForm`](@ref):
`error_control_variables = water_content` evaluates `water_content` on pressure head
rather than selecting the stored water-content block. Use `nothing` to retain error
control in the evolved variables.

The converted residual is reduced with `Trixi.ode_norm` (an MPI-aware RMS norm).
`internalnorm` remains available to OrdinaryDiffEq for its original estimator and other
solver operations; it is not applied to the converted residual. A controller wrapper
replaces the original error estimate before step-size selection and acceptance. The
original error scaling and norm are still computed;
stages, Jacobians, and linear solves are not repeated. The option has no effect when
`adaptive = false`.

Converted error control supports the in-place implementations of `Rodas4`, `Rodas42`,
`Rodas4P`, `Rodas4P2`, `Rodas5`, `Rodas5P`, `Rodas5Pe`, and `Rodas6P` with the default
step limiter. It wraps an OrdinaryDiffEq controller, such as
`OrdinaryDiffEqCore.PIController(algorithm)`; if `controller` is `nothing`, it creates
that PI controller. Unsupported algorithms are rejected when conversion is enabled.

Without conversion, HydroTrixi.jl's default controller is used only with `Rodas5P`.
For any other integration algorithm, OrdinaryDiffEq.jl selects its default controller
unless `controller` is passed explicitly.
"""
function solve_implicit(ode::SciMLBase.ODEProblem{U, T, I, P},
                        algorithm = default_algorithm(ode);
                        dt,
                        adaptive = true,
                        abstol = 1.0e-11,
                        reltol = 1.0e-7,
                        controller = default_stepsize_controller(algorithm, ode),
                        dtmin = zero(last(ode.tspan) - first(ode.tspan)),
                        dtmax = last(ode.tspan) - first(ode.tspan),
                        force_dtmin = false,
                        failfactor = 2, # not used by Rodas5P (a linearly implicit method)
                        maxiters = typemax(Int),
                        internalnorm = evolved_variable_norm(ode.p),
                        error_control_variables = nothing,
                        save_everystep = false,
                        unstable_check = Trixi.mpi_isparallel() ?
                                         Trixi.ode_unstable_check :
                                         DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK,
                        kwargs...) where {U, T, I, P <: SemidiscretizationImplicit}
    if adaptive && error_control_variables !== nothing
        if error_control_uses_evolved_variables(error_control_variables,
                                                ode.p.operator_temporal)
            internalnorm = evolved_variable_norm(ode.p)
        else
            controller = StepsizeControllerMappedError(controller, algorithm, ode,
                                                        error_control_variables; reltol)
        end
    end
    common_options = (; dt, adaptive, dtmin, dtmax, force_dtmin, failfactor,
                      maxiters, internalnorm, save_everystep, unstable_check)
    adaptive_options = adaptive ? (; abstol, reltol, controller) : (;)
    options = merge(common_options, adaptive_options, (; kwargs...))
    return SciMLBase.solve(ode, algorithm; options...)
end
