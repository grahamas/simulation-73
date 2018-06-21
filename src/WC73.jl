module WC73

using Parameters
using Simulation

@with_kw struct WC73{T<:Real} <: Model
    α::Array{T}
    β::Array{T}
    τ::Array{T}
    space::Space{T}
    connectivity::Array{Connectivity{T}}
    nonlinearity::Array{Nonlinearity{T}}
    stimulus::Array{Stimulus{T}}
    WC73{T<:Real}(α, β, τ, s, c, n, t) = new(α, β, τ, s, c, n, t)
end

# * Calculated WC73 Simulation Type
mutable struct CalculatedWC73{T,C<:Connectivity,
                      N<:Nonlinearity,S<:Stimulus}
    α::Array{T}
    β::Array{T}
    τ::Array{T}
    connectivity::Array{Calculated{C{T}}}
    nonlinearity::Array{Calculated{N{T}}}
    stimulus::Array{Calculated{S{T}}}
end

function CalculatedWC73(wc::WC73{T,C,N,S}) where {T<:Real, C<:Connectivity,
                                                  N<:Nonlinearity, S<:Stimulus}
    CalculatedWC73{T,C,N,S}(
        wc.α, wc.β, wc.τ,
        Calculated.(wc.connectivity, wc.space),
        Calculated.(wc.nonlinearity),
        Calculated.(wc.stimulus,wc.space))
end

function update_from_p!(cwc::CalculatedWC73{C,N,S}, new_p, p_search::ParameterSearch{WC73})
    new_model = model_from_p(p_search, new_p)
    cwc.α = new_model.α
    cwc.β = new_model.β
    cwc.τ = new_model.τ
    update!(cwc, new_model.connectivity)
    update!(cwc, new_model.nonlinearity)
    update!(cwc, new_model.stimulus)
end
update!(cwc::CalculatedWilsonCowan73, c::Connectivity) = update!(cwc.connectivity, c)
update!(cwc::CalculatedWilsonCowan73, n::Nonlinearity) = update!(cwc.nonlinearity, n)
update!(cwc::CalculatedWilsonCowan73, s::Stimulus) = update!(cwc.stimulus, s)

function get_values(cwc::CalculatedWilsonCowan73)
    (cwc.α, cwc.β, cwc.τ, cwc.connectivity.value, cwc.nonlinearity.value, cwc.stimulus.value)
end

# * Problem generation
function make_problem_generator(p_search::ParameterSearch{WC73})
    model = initial_model(p_search)
    tspan = time_span(p_search)

    u0 = initial_value(model)

    n_pops = length(model.α)
    cwc = CalculatedWilsonCowan73(model)
    function problem_generator(prob, new_p)
        update_from_p!(cwc, new_p, p_search)
        α, β, τ, connectivity_fn, nonlinearity_fn, stimulus_fn = get_values(cwc)
        function WilsonCowan73!(dA, A, p, t)
            for i in 1:n_pops
                dA[i,:] .= (-α[i] .* A[i,:] .+ β[i] .* (1.-A[i]) .* nonlinearity_fn[i](sum(connectivity_fn[i,j](A[j,:]) for j in 1:n_pops) .+ stimulus_fn[i](t))) ./ τ[i]
            end
        end
        ODEProblem(WilsonCowan73!, u0, tspan, new_p)
    end
    initial_p = param_vector(cwc)
    initial_problem = problem_generator(nothing, initial_p)
    return initial_problem, problem_generator
end



end
