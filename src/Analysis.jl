module Analysis
using CalculatedParameters, Parameters
using Parameters
using DiffEqBase: AbstractODESolution, ODESolution
using OrdinaryDiffEq: ODECompositeSolution, InterpolationData, CompositeInterpolationData
using Modeling
using Records

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots


# * Analysis types
# TODO: DANGER AbstractFigure is implemented in Recipes(Base?)
abstract type AbstractResults{T,N,M<:Model} end
abstract type AbstractPlotSpecification end
function output_name(af::AbstractPlotSpecification)
	af.output_name
end

function plot_and_save(plot_spec::APS, results::AbstractResults, output::AbstractOutput) where {APS <: AbstractPlotSpecification}
	save_fn(fn, plt) = savefig(plt, fn)
	output(save_fn, output_name(plot_spec), plot(plot_spec, results; plot_spec.kwargs...))
end

struct SubSampler{T}
	dt::T
	space_strides::Array{Int,1}
	time_window::Tuple{T,T}
	space_window::Tuple{T,T}
end
function SubSampler(; dt::T, space_strides::Array{Int}) where T
    SubSampler{T}(dt, space_strides, (0.0, Inf), (-Inf, Inf))
end

tail_args(_, rest...) = rest
tail(tup::Tuple) = tail_args(tup...)
tail(arr::Array{T,1}) where T = tail_args(arr...)
subinds(steps::AbstractArray{T,1}, inds::Tuple{R}) where {T,S,R <: AbstractArray{S}} = (first(inds)[1:first(steps):end], subinds(tail(steps), tail(inds))...)
subinds(::Tuple{}, ::Tuple{}) = ()
subinds(::Tuple{}, ::Any) = error("Too many spatial strides in subsampling.")
subinds(::Any, ::Tuple{}) = error("Too few spatial strides in subsampling.")

function sample end

@with_kw struct Analyses{T}
	subsampler::Union{SubSampler{T},Nothing} = nothing
	plots::Array{AbstractPlotSpecification}
end

@with_kw struct Results{T,N,M,uType} <: AbstractResults{T,N,M}
	model::M
	solution::ODES where {F,uType2,tType<:Array{T,1},DType,kType,cacheType,rateType,P,A,IType <: Union{InterpolationData{F,uType2,tType,kType,cacheType},CompositeInterpolationData{F,uType,tType,kType,cacheType}}, ODES <: Union{ODESolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType},ODECompositeSolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType}}}
end
# @with_kw struct SubSampledResults{T,N,M} <: AbstractResults{T,N,M}
# 	model::M
# 	solution::DESolution
# 	subsampler::SubSampler
# end
@with_kw struct SubSampledResults{T,N,M,uType} <: AbstractResults{T,N,M}
	model::M
	solution::ODES where {F,uType2,tType<:Array{T,1},DType,kType,cacheType,rateType,P,A,IType <: Union{InterpolationData{F,uType2,tType,kType,cacheType},CompositeInterpolationData{F,uType,tType,kType,cacheType}}, ODES <: Union{ODESolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType},ODECompositeSolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType}}}
	subsampler::SubSampler{T}
end

function SubSampledResults(model::M, solution::ODECompositeSolution{T,N,a,b,c,d,e,f,g,IType}, subsampler::SubSampler{T}) where {T,N,M,IType,a,b,c,d,e,f,g}
	SubSampledResults{T,N,M,IType}(model, solution, subsampler)
end

function Results(model::M, solution::AbstractODESolution{T,N}, subsampler::Nothing) where {T,N,M <: Model{T}}
	Results{T,N,M}(model, solution)
end
function Results(model::M, solution::AbstractODESolution{T,N}, subsampler::SubSampler) where {T,N,M <: Model{T}}
	SubSampledResults(model, solution, subsampler)
end

function resample(r::AbstractResults{T,N,M}, s::SubSampler) where {T,N,M<:Model{T}}
	Results(r.model, r.solution, s)
end
function resample(r::SubSampledResults{T,N,M}; dt=nothing, space_strides=nothing,
		time_window=nothing, space_window=nothing) where {T,N,M<:Model{T}}
	if time_window == nothing
		time_window = r.subsampler.time_window
	end
	if space_window == nothing
		space_window = r.subsampler.space_window
	end
	if dt == nothing
		dt = r.subsampler.dt
	end
	if space_strides == nothing
		space_strides = r.subsampler.space_strides
	end
	resample(r, SubSampler(dt, space_strides, time_window, space_window))
end

function get_space_dx(results::Results)
	Colon()
end
function get_space_dx(results::SubSampledResults)
	frame = get_space(results.model)
	windowed_indices = axes(frame)[1][results.subsampler.space_window[1] .<= frame .<= results.subsampler.space_window[end]]
	subinds(results.subsampler.space_strides, (windowed_indices,))[1] # ASSUMES D = 1!!!
end

function sample_space(frame::AbstractArray{T}, results::SubSampledResults) where {T,M}
	getindex(frame, get_space_dx(results), Colon()) # Assumes 1 trailing pop D
end
function sample_space(frame, results::Results)
	frame
end

function Modeling.get_space(results::AbstractResults)
	sample_space(get_space(results.model), results)
end

function get_time(results::Results)
	results.solution.t
end
function get_time(results::SubSampledResults)
	0:results.dt:results.solution.t[end]
end

function get_pop(results::AbstractResults, pop_num::Int)
	t_dx = get_time(results)
	x_dx = get_space(results)
	results.solution(t_dx, idxs=CartesianIndex.(get_space_dx(results), pop_num))
end

Base.iterate(r::Results, state...) = iterate(tuples(r.solution), state...)
function Base.iterate(r::SR) where {T, N, M <: Model{T}, SR <: SubSampledResults{T,N,M}}
	start = r.subsampler.time_window[1]
	sampled = sample_space(r.solution(start), r)
	((sampled, start), start)
end
function Base.iterate(r::SR, prev_t) where {T, N, M <: Model{T}, SR <: SubSampledResults{T,N,M}}
	new_t = prev_t + r.subsampler.dt
	if new_t > r.solution.t[end] || new_t > r.subsampler.time_window[end]
		return nothing
	else
		sampled = sample_space(r.solution(new_t), r)
		return ((sampled, new_t), new_t)
	end
end

function analyse(a::Analyses, results::AbstractResults, output::AbstractOutput)
    @info "Begin analysis."
    #return results
    a.plots .|> (plot_type) -> plot_and_save(plot_type, results, output)
end

export AbstractResults, AbstractPlotSpecification, output_name
export plot_and_save
export SubSampler, sample, Analyses, Results, SubSampledResults
export spatiotemporal_data
export analyse
export get_space_dx, get_space, get_time, resample, get_pop

end
