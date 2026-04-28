@doc raw"""
    SemidiscretizationConstitutive{Semidiscretization, ConstitutiveRelation}

A semidiscretization wrapper that augments an internal semidiscretization with a
constitutive relation ``u_\mathrm{evolved} = \Theta(u_\mathrm{state})`` to form a mixed system of the form
```math
\partial_t u_\mathrm{evolved} = R(u_\mathrm{state}, t), \qquad
0 = u_\mathrm{evolved} - \Theta(u_\mathrm{state}).
```
"""
struct SemidiscretizationConstitutive{Semidiscretization, ConstitutiveRelation} <:
       Trixi.AbstractSemidiscretization
    semidiscretization::Semidiscretization
    constitutive_relation::ConstitutiveRelation
end

@inline Base.ndims(semi::SemidiscretizationConstitutive) = ndims(semi.semidiscretization)
@inline Base.real(semi::SemidiscretizationConstitutive) = real(semi.semidiscretization)
@inline Trixi.nvariables(semi::SemidiscretizationConstitutive) = 2 *
                                                                 Trixi.nvariables(semi.semidiscretization)

@inline function _check_even_length(u_ode)
    iseven(length(u_ode)) ||
        throw(ArgumentError("Expected an even number of constitutive degrees of freedom."))
    return nothing
end

# Create views into the state and evolved variable blocks of the ODE vector
@inline function split_variable_views(u_ode)
    half = length(u_ode) ÷ 2
    return @view(u_ode[1:half]), @view(u_ode[(half + 1):end])
end

@inline evolved_variable_view(u_ode) = @view(u_ode[1:(length(u_ode) ÷ 2)])
@inline state_variable_view(u_ode) = @view(u_ode[(length(u_ode) ÷ 2 + 1):end])

function Trixi.compute_coefficients(t, semi::SemidiscretizationConstitutive)
    # First compute the state variable coefficients from the internal semidiscretization
    state_coefficients = Trixi.compute_coefficients(t, semi.semidiscretization)
    evolved_coefficients = similar(state_coefficients)
    equations = semi.semidiscretization.equations

    # Apply the constitutive relation to compute the evolved variable coefficients
    @inbounds for i in eachindex(evolved_coefficients, state_coefficients)
        evolved_coefficients[i] = semi.constitutive_relation(state_coefficients[i],
                                                             equations)
    end

    # Concatenate into a single vector
    return vcat(evolved_coefficients, state_coefficients)
end

function rhs_constitutive!(du_ode, u_ode, semi::SemidiscretizationConstitutive, t)
    evolved_variable, state_variable = split_variable_views(u_ode)
    evolved_variable_rhs, state_variable_rhs = split_variable_views(du_ode)

    # First block of du_ode: R(u_state, t)
    Trixi.default_rhs(semi.semidiscretization)(evolved_variable_rhs, state_variable,
                                               semi.semidiscretization, t)
    (; equations) = semi.semidiscretization

    # Second block of du_ode: u_evolved - constitutive_relation(u_state)
    @inbounds for i in eachindex(state_variable_rhs, evolved_variable, state_variable)
        state_variable_rhs[i] = evolved_variable[i] -
                                semi.constitutive_relation(state_variable[i], equations)
    end

    return nothing
end

function Trixi.semidiscretize(semi::SemidiscretizationConstitutive, tspan;
                              reset_threads = true)
    if reset_threads
        Trixi.Polyester.reset_threads!()
    end

    u0_ode = Trixi.compute_coefficients(first(tspan), semi)
    _check_even_length(u0_ode)

    # Set mass matrix to [I, 0; 0, 0]
    half = length(u0_ode) ÷ 2
    diagonal_entries = zeros(eltype(u0_ode), length(u0_ode))
    @inbounds diagonal_entries[1:half] .= one(eltype(u0_ode))
    mass_matrix = Diagonal(diagonal_entries)

    ode_function = SciMLBase.ODEFunction{true, SciMLBase.FullSpecialize}(rhs_constitutive!;
                                                                         mass_matrix = mass_matrix)
    return SciMLBase.ODEProblem{true, SciMLBase.FullSpecialize}(ode_function,
                                                                u0_ode,
                                                                tspan,
                                                                semi)
end
