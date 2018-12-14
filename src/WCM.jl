module WCM

using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update!
using Simulating
using Modeling
using WCMConnectivity
using WCMNonlinearity
using WCMStimulus
using Meshes
using Exploration
using DifferentialEquations
using Targets
using Records
import Records: required_modules

struct WCMSpatial1D{T,C<:Connectivity{T},
                            L<:Nonlinearity{T},S<:Stimulus{T},SP<:Space{T}} <: Model{T}
    α::Array{T,1}
    β::Array{T,1}
    τ::Array{T,1}
    space::SP
    connectivity::Array{C,2}
    nonlinearity::Array{L,1}
    stimulus::Array{S,1}
    pop_names::Array{<:AbstractString}
end

function WCMSpatial1D(; pop_names::Array{<:AbstractString,1}, α::Array{T,1}, β::Array{T,1}, τ::Array{T,1},
        space::SP, connectivity::Array{C,2}, nonlinearity::Array{L,1}, stimulus::Array{S,1}) where {T,C<:Connectivity{T},L<:Nonlinearity{T},S<:Stimulus{T},SP<:Space{T}}
    WCMSpatial1D{T,C,L,S,SP}(α,β,τ,space,connectivity,nonlinearity,stimulus,pop_names)
end


space_array(model::WCMSpatial1D) = Calculated(model.space).value

export WCMSpatial1D, space_array

import Exploration: base_type
function base_type(::Type{WCMSpatial1D{T1,T2,T3,T4}}) where {T1,T2,T3,T4}
    BT1 = base_type(T1); BT2 = base_type(T2); BT3 = base_type(T3)
    BT4 = base_type(T4)
    return WCMSpatial1D{BT1,BT2,BT3,BT4}
end
export base_type

# * Calculated WC73 Simulation Type
struct CalculatedWCMSpatial1D{T,C,L,S,CC<:CalculatedParam{C},CL <: CalculatedParam{L},CS <: CalculatedParam{S}} <: CalculatedParam{WCMSpatial1D{T,C,L,S}}
    α::Array{T,1}
    β::Array{T,1}
    τ::Array{T,1}
    connectivity::Array{CC,2}
    nonlinearity::Array{CL,1}
    stimulus::Array{CS,1}
    function CalculatedWCMSpatial1D{T,C,L,S,CC,CL,CS}(α::Array{T},
                                                β::Array{T},
                                                τ::Array{T},
                                                connectivity::Array{CC},
                                                nonlinearity::Array{CL},
                                                stimulus::Array{CS}) where {T,C,L,S,CC,CL,CS}
        new(α,β,τ,connectivity,nonlinearity,stimulus)
    end
end

function CalculatedWCMSpatial1D(wc::WCMSpatial1D{T,C,L,S}) where {T<:Real,
                                                  C<:Connectivity{T},
                                                  L<:Nonlinearity{T},
                                                  S<:Stimulus{T}}
    connectivity = Calculated.(wc.connectivity, Ref(wc.space))
    nonlinearity = Calculated.(wc.nonlinearity)
    stimulus = Calculated.(wc.stimulus,Ref(wc.space))
    CC = eltype(connectivity)
    CL = eltype(nonlinearity)
    CS = eltype(stimulus)
    CalculatedWCMSpatial1D{T,C,L,S,CC,CL,CS}(
        wc.α, wc.β, wc.τ,
        connectivity, nonlinearity, stimulus)
end

function Calculated(α::Array{T}, β::Array{T}, τ::Array{T},
    connectivity::Array{CC}, nonlinearity::Array{CL}, stimulus::Array{CS}) where {T,C<:Connectivity{T},L<:Nonlinearity{T},S<:Stimulus{T},CC<:CalculatedParam{C},CL<:CalculatedParam{L},CS<:CalculatedParam{S}}
    CalculatedWCMSpatial1D{T,C,L,S,CC,CL,CS}(α, β, τ, connectivity, nonlinearity, stimulus)
end

function Calculated(wc::WCMSpatial1D)
    CalculatedWCMSpatial1D(wc)
end

function update_from_p!(cwc::CalculatedWCMSpatial1D{<:Real}, new_p, p_search::ParameterSearch{<:WCMSpatial1D})
    # Use the variable model stored by p_search to create static model
    new_model = model_from_p(p_search, new_p)

    # Update the calculated values from the new static model
    cwc.α = new_model.α
    cwc.β = new_model.β
    cwc.τ = new_model.τ
    update!(cwc.connectivity, new_model.connectivity, space)
    update!(cwc.nonlinearity, new_model.nonlinearity)
    update!(cwc.stimulus, new_model.stimulus, space)
end

function get_values(cwc::CalculatedWCMSpatial1D)
    (cwc.α, cwc.β, cwc.τ, get_value.(cwc.connectivity), get_value.(cwc.nonlinearity), get_value.(cwc.stimulus))
end

import Exploration: make_problem_generator
# * Problem generation
function make_problem_generator(p_search::ParameterSearch{<:WCMSpatial1D})
    model = initial_model(p_search)
    tspan = time_span(p_search)

    u0 = initial_value(model)

    n_pops = length(model.α)
    cwc = CalculatedWCMSpatial1D(model)
    function problem_generator(prob, new_p)
        update_from_p!(cwc, new_p, p_search)
        α, β, τ, P, connectivity_mx, nonlinearity_fn, stimulus_fn = get_values(cwc)
        function WilsonCowan73!(dA::Array{T,2}, A::Array{T,2}, p::Array{T,1}, t::T)::Nothing where {T<:Float64}
            for i in 1:n_pops
                stim_val::Array{T,1} = stimulus_fn[i](t)
                nonl_val::Array{T,1} = nonlinearity_fn[i](sum(connectivity_mx[i,j]::Array{T,2} * A[:,j] for j in 1:n_pops) .+ stim_val)
                dA[:,i] .= (-α[i] .* A[:,i] .+ β[i] .* (1.0 .- A[:,i]) .*  nonl_val + P[i]) ./ τ[i]
            end
        end
        ODEProblem(WilsonCowan73!, u0, tspan, new_p)
    end
    initial_problem = problem_generator(nothing, p_search.initial_p)
    return initial_problem, problem_generator
end
export make_problem_generator

function generate_problem(simulation::Simulation{T,<:WCMSpatial1D}) where T
    tspan = time_span(simulation)
    model = simulation.model
    u0 = initial_value(model)
    n_pops = length(model.pop_names)
    cwc = CalculatedWCMSpatial1D(model)
    α, β, τ, connectivity_mx, nonlinearity_objs, stimulus_objs = get_values(cwc)
    function WilsonCowan73!(dA::Array{T,2}, A::Array{T,2}, p::Union{Array{T,1},Nothing}, t::T)::Nothing where {T<:Float64}
        for i in 1:n_pops
            stim_val::Array{T,1} = stimulus(stimulus_objs[i], t) # I'll bet it goes faster if we pull this out of the loop
            nonl_val::Array{T,1} = nonlinearity.(Ref(nonlinearity_objs[i]), sum(connectivity_mx[i,j]::Array{T,2} * A[:,j] for j in 1:n_pops) .+ stim_val)
            dA[:,i] .= (-α[i] .* A[:,i] .+ β[i] .* (1.0 .- A[:,i]) .*  nonl_val) ./ τ[i]
        end
    end
    return ODEProblem(WilsonCowan73!, u0, tspan, nothing)
end

export generate_problem

end
