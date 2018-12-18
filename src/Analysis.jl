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
abstract type AbstractResults{T,N,M<:Model{T,N},uType} end # TODO specify uType in terms of T and N
abstract type AbstractPlotSpecification end
function output_name(af::AbstractPlotSpecification)
	af.output_name
end

function plot_and_save(plot_spec::APS, results::AbstractResults, output::AbstractOutput) where {APS <: AbstractPlotSpecification}
	@info "Entered plot and save"
	save_fn(fn, plt) = savefig(plt, fn)
	@info "Calling output plot"
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

struct Results{T,N,M,uType,O} <: AbstractResults{T,N,M,uType}
	model::M
	solution::O# where {F,uType2,tType<:Array{T,1},DType,kType,cacheType,rateType,P,A,IType <: Union{InterpolationData{F,uType,tType,kType,cacheType},CompositeInterpolationData{F,uType,tType,kType,cacheType}}, O <: Union{ODESolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType},ODECompositeSolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType}}}
end

struct SubSampledResults{T,N,M,uType,ODES} <: AbstractResults{T,N,M,uType}
	model::M
	solution::ODES
	subsampler::SubSampler{T}
end

function SubSampledResults(model::M, solution::O, subsampler::SubSampler{T}) where {T,N,NSol,M<:Model{T,N},uType,F,uType2,tType<:Array{T,1},DType,kType,cacheType,rateType,P,A,IType <: Union{InterpolationData{F,uType,tType,kType,cacheType},CompositeInterpolationData{F,uType,tType,kType,cacheType}}, O <: Union{ODESolution{T,NSol,uType,uType2,DType,tType,rateType,P,A,IType},ODECompositeSolution{T,NSol,uType,uType2,DType,tType,rateType,P,A,IType}}}
	SubSampledResults{T,N,M,uType,O}(model, solution, subsampler)
end
function Results(model::M, solution::AbstractODESolution{T,N}, subsampler::SubSampler) where {T,N,M <: Model{T}}
	@info "General results constructor"
	SubSampledResults(model, solution, subsampler)
end

function resample(r::AbstractResults{T,N,M,uType}, s::SubSampler{T}) where {T,N,M<:Model{T},uType}
	@info "Resampling abstract with subsampler"
	Results(r.model, r.solution, s)
end
function resample(r::SubSampledResults{T,N,M,uType}; dt=nothing, space_strides=nothing,
		time_window=nothing, space_window=nothing) where {T,N,M<:Model{T,N},uType}
	@info "Resampling kw"
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
	@info "Going into next resampling"
	resample(r, SubSampler(dt, space_strides, time_window, space_window))
end

function get_space_dx(results::Results)
	Colon()
end
function get_space_dx(results::SubSampledResults{T,N,M,uType}) where {T,N,M,uType}
	frame = get_space(results.model)
	windowed_indices = axes(frame)[1][results.subsampler.space_window[1] .<= frame .<= results.subsampler.space_window[end]]
	subinds(results.subsampler.space_strides, (windowed_indices,))[1] # ASSUMES D = 1!!!
end

function sample_space(frame::AbstractArray{T}, results::SubSampledResults{T,N,M,uType,ODES})::frameType where {T,N,M,frameType,uType<:Array{frameType,1},ODES}
	getindex(frame, get_space_dx(results), Colon()) # Assumes 1 trailing pop D
end
function sample_space(frame, results::Results)
	frame
end

function Modeling.get_space(results::AbstractResults{T,N,M,uType}) where {T,N,M,uType}
	sample_space(get_space(results.model), results)
end

function get_time(results::Results{T}) where T
	results.solution.t
end
function get_time(results::SubSampledResults{T}) where T
	0:results.dt:results.solution.t[end]
end

function get_pop(results::AbstractResults, pop_num::Int)
	t_dx = get_time(results)
	x_dx = get_space(results)
	results.solution(t_dx, idxs=CartesianIndex.(get_space_dx(results), pop_num))
end

Base.iterate(r::Results{T,N,M,uType}, state...) where {T,N,M,uType} = iterate(tuples(r.solution), state...)
function Base.iterate(r::SR) where {T, N, M <: Model{T,N}, uType, SR <: SubSampledResults{T,N,M,uType}}
	start = r.subsampler.time_window[1]
	sampled = sample_space(r.solution(start), r)
	((sampled, start), start)
end
function Base.iterate(r::SR, prev_t) where {T, N, M <: Model{T,N}, uType, SR <: SubSampledResults{T,N,M,uType}}
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
    a.plots .|> (plot_spec) -> plot_and_save(plot_spec, results, output)
end

export AbstractResults, AbstractPlotSpecification, output_name
export plot_and_save
export SubSampler, sample, Analyses, Results, SubSampledResults
export spatiotemporal_data
export analyse
export get_space_dx, get_space, get_time, resample, get_pop

end
