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
    return OrdinaryDiffEqRosenbrock.Rodas5P(; autodiff, standardtag, concrete_jac,
                                            diff_type, linsolve, precs,
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
    return OrdinaryDiffEqCore.NewPIController(controller_type, algorithm;
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

@doc raw"""
    solve_implicit(ode, algorithm=default_algorithm(ode);
                   dt, adaptive=true, abstol=1.0e-11, reltol=1.0e-7,
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
between `1.0` and `1.2`. These controller parameters reproduce the defaults in
`OrdinaryDiffEqCore` v3.33.1. For a time integration method of order ``p=5``, the
[order-dependent defaults](https://github.com/SciML/OrdinaryDiffEq.jl/blob/3eb62b46769c5db70c131f6b4331ebb0a6864117/lib/OrdinaryDiffEqCore/src/alg_utils.jl#L526-L550)
give ``\beta_1=7/(10p)=0.14`` and ``\beta_2=2/(5p)=0.08``, as well as the safety factor
and steady-step deadband used here. The clipping factors come from the
[step-size-factor defaults](https://github.com/SciML/OrdinaryDiffEq.jl/blob/3eb62b46769c5db70c131f6b4331ebb0a6864117/lib/OrdinaryDiffEqCore/src/alg_utils.jl#L232-L242),
whereas the first-step growth factor and stored previous error come from the
[`NewPIController` constructor](https://github.com/SciML/OrdinaryDiffEq.jl/blob/3eb62b46769c5db70c131f6b4331ebb0a6864117/lib/OrdinaryDiffEqCore/src/integrators/controllers.jl#L384-L398).

The defaults also use [`state_variable_norm`](@ref), disable saving every accepted step,
and allow `typemax(Int)` iterations.

HydroTrixi.jl's default controller is used only with `Rodas5P`. For any other integration
algorithm, OrdinaryDiffEq.jl selects its default controller unless `controller` is passed
explicitly.
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
                        internalnorm = state_variable_norm(ode.p),
                        save_everystep = false,
                        unstable_check = Trixi.mpi_isparallel() ?
                                         Trixi.ode_unstable_check :
                                         DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK,
                        kwargs...) where {U, T, I, P <: SemidiscretizationImplicit}
    common_options = (; dt, adaptive, dtmin, dtmax, force_dtmin, failfactor, maxiters,
                      internalnorm, save_everystep, unstable_check)
    adaptive_options = adaptive ? (; abstol, reltol, controller) : (;)
    options = merge(common_options, adaptive_options, (; kwargs...))
    return SciMLBase.solve(ode, algorithm; options...)
end
