module WC73

using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update!
using Modeling
using Analysis
using WCMAnalysis
using WCMConnectivity
using WCMNonlinearity
using WCMStimulus
using Meshes
using Exploration
using DifferentialEquations
using Targets
using Records
import Records: required_modules

@with_kw struct WCMSpatial1D{T,C<:Connectivity{T},
                            L<:Nonlinearity{T},S<:Stimulus{T}} <: Model{T}
    α::Array{T}
    β::Array{T}
    τ::Array{T}
    space::Space{T}
    connectivity::Array{C}
    nonlinearity::Array{L}
    stimulus::Array{S}
    function WCMSpatial1D{T,C,L,S}(α::Array{T},
        β::Array{T}, τ::Array{T},
        s::Space{T}, c::Array{C},
        n::Array{L}, t::Array{S}) where {T,
                                         C<:Connectivity{T},
                                         L<:Nonlinearity{T},
                                         S<:Stimulus{T}}
        new(α, β, τ, s, c, n, t)
    end
end

function required_modules(::Type{M}) where {M <: WCMSpatial1D}
    [Modeling, WC73, Meshes, CalculatedParameters, WCMConnectivity, WCMNonlinearity, WCMStimulus]
end

export WCMSpatial1D, required_modules

import Exploration: base_type
function base_type(::Type{WCMSpatial1D{T1,T2,T3,T4}}) where {T1,T2,T3,T4}
    BT1 = base_type(T1); BT2 = base_type(T2); BT3 = base_type(T3)
    BT4 = base_type(T4)
    return WCMSpatial1D{BT1,BT2,BT3,BT4}
end
export base_type

# * Calculated WC73 Simulation Type
mutable struct CalculatedWCMSpatial1D{T,C,L,S, CC<:CalculatedParam{C},CL <: CalculatedParam{L},CS <: CalculatedParam{S}} <: CalculatedParam{WCMSpatial1D{T,C,L,S}}
    α::Array{T}
    β::Array{T}
    τ::Array{T}
    connectivity::Array{CC}
    nonlinearity::Array{CL}
    stimulus::Array{CS}
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
    connectivity = Calculated.(wc.connectivity, wc.space)
    nonlinearity = Calculated.(wc.nonlinearity)
    stimulus = Calculated.(wc.stimulus,wc.space)
    CC = eltype(connectivity)
    CL = eltype(nonlinearity)
    CS = eltype(stimulus)
    CalculatedWCMSpatial1D{T,C,L,S,CC,CL,CS}(
        wc.α, wc.β, wc.τ,
        connectivity, nonlinearity, stimulus)
end

function update_from_p!(cwc::CalculatedWCMSpatial1D{<:Real}, new_p, p_search::ParameterSearch{<:WCMSpatial1D})
    new_model = model_from_p(p_search, new_p)
    cwc.α = new_model.α
    cwc.β = new_model.β
    cwc.τ = new_model.τ
    update!(cwc, new_model.connectivity)
    update!(cwc, new_model.nonlinearity)
    update!(cwc, new_model.stimulus)
end
update!(cwc::CalculatedWCMSpatial1D, c::Array{<:Connectivity}) = update!(cwc.connectivity, c)
update!(cwc::CalculatedWCMSpatial1D, n::Array{<:Nonlinearity}) = update!(cwc.nonlinearity, n)
update!(cwc::CalculatedWCMSpatial1D, s::Array{<:Stimulus}) = update!(cwc.stimulus, s)

function get_values(cwc::CalculatedWCMSpatial1D)
    get_value = (el) -> el.value
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
        α, β, τ, connectivity_mx, nonlinearity_fn, stimulus_fn = get_values(cwc)
        function WilsonCowan73!(dA::Array{T,2}, A::Array{T,2}, p::Array{T,1}, t::T) where {T<:Float64}
            for i in 1:n_pops
                #stim_val::Array{T,1} = stimulus_fn[i](t)
                #nonl_val::Array{T,1} = nonlinearity_fn[i](sum(connectivity_mx[i,j] * A[:,j] for j in 1:n_pops) .+ stim_val)
                dA[:,i] .= (-α[i] .* A[:,i] .+ β[i] .* (1.-A[:,i]) .*  nonlinearity_fn[i](sum(connectivity_mx[i,j] * A[:,j] for j in 1:n_pops) .+ stimulus_fn[i](t))) ./ τ[i]
            end
        end
        ODEProblem(WilsonCowan73!, u0, tspan, new_p)
    end
    initial_problem = problem_generator(nothing, p_search.initial_p)
    return initial_problem, problem_generator
end
export make_problem_generator

function analyse(soln::DESolution, write_fn::Function, model::WCMSpatial1D{T};
                 down_sampling=nothing, nonlinearity=nothing,
                 pop_names=nothing, activity_gif=nothing, heatmap=nothing) where T <: Real
    space = calculate(model.space)
    timeseries = soln.u
    n_pops = length(model.α)
    n_space = length(space)
    @assert all(size(timeseries, [1,2]) .== (n_space, n_pops))
    if (down_sampling != nothing) ds_t, ds_x, ds_timeseries = down_sample(soln.t, space, timeseries; down_sampling...) end
    if (heatmap != nothing) plot_heatmap(ds_t, ds_x, ds_timeseries, write_fn; heatmap...) end
    #pop_peak_timeseries = calc_pop_peak_timeseries(timeseries, 0)
    #if (nonlinearity != nothing) plot_nonlinearity(soln.prob.p.nonlinearity_fn, write_fn, pop_names; nonlinearity) end
    if (activity_gif != nothing) plot_activity_gif(ds_t, ds_x, ds_timeseries, write_fn; activity_gif...) end
end

export analyse



function Analysis.analyse(soln::DESolution, output::Output, model::WCMSpatial1D; sampling=nothing,
                          nonlinearity=nothing, pop_names=nothing, activity_gif=nothing, heatmap=nothing)
    ds_t, ds_x, ds_timeseries = sample_timeseries(soln, model; sampling...)
    if (heatmap != nothing) plot_heatmap(ds_t, ds_x, ds_timeseries, output; heatmap...) end
    #pop_peak_timeseries = calc_pop_peak_timeseries(timeseries, 0)
    #if (nonlinearity != nothing) plot_nonlinearity(soln.prob.p.nonlinearity_fn, output, pop_names; nonlinearity) end
    if (activity_gif != nothing) plot_activity_gif(ds_t, ds_x, ds_timeseries, output; activity_gif...) end
end

function Analysis.analyse(sim::Simulation{WCMSpatial1D}, soln)
    analyse(soln, sim.output, sim.model; sim.analyses...)
end

export analyse

end
