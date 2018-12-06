module Analysis
using CalculatedParameters, Parameters
using Parameters
using DiffEqBase
using Modeling
using Records

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots


# * Analysis types
# TODO: DANGER AbstractFigure is implemented in Recipes(Base?)
abstract type AbstractResults{M<:Model} end
abstract type AbstractPlotSpecification end
function output_name(af::AbstractPlotSpecification)
	af.output_name
end

function plot_and_save(plot_spec::APS, results::AbstractResults, output::AbstractOutput) where {APS <: AbstractPlotSpecification}
	save_fn(fn, plt) = savefig(plt, fn)
	output(save_fn, output_name(plot_spec), plot(plot_spec, results; plot_spec.kwargs...))
end



@with_kw struct SubSampler
	dt::Float64
	space_strides::Array{Integer}
end
function SubSampler(dt::Float64, spatial_stride::Int)
    SubSampler(dt, [spatial_stride])
end

subinds(steps, inds) = (first(inds)[1:first(steps):end], subinds(tail(steps), tail(inds))...)
subinds(::Tuple{}, ::Tuple{}) = ()
subinds(::Tuple{}, ::Any) = error("Too many spatial strides in subsampling.")
subinds(::Any, ::Tuple{}) = error("Too few spatial strides in subsampling.")

function sample end

@with_kw struct Analyses
	subsampler::Union{SubSampler,Nothing} = nothing
	plots::Array{AbstractPlotSpecification}
end

@with_kw struct Results{M} <: AbstractResults{M}
	model::M
	solution::DESolution
end
@with_kw struct SubSampledResults{M} <: AbstractResults{M}
	model::M
	solution::DESolution
	subsampler::SubSampler
end

function get_space_dx(results::Results)
	Colon()
end
function get_space_dx(results::SubSampledResults)
	frame = results.solution.x
	subinds(results.subsampler.space_strides, indices(frame))
end
function sample_space(frame::Array{T}, results::SubSampledResults) where {T,M}
	getindex(frame, get_space_dx(results))
end

function get_space(results::Results)
	results.solution.x
end
function get_space(results::SubSampledResults)
	sample_space(results.solution.x, results)
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

export get_space_dx, get_space, get_time, get_pop

function Results(model::M, solution::DESolution, subsampler::Nothing) where {M <: Model}
	Results(model, solution)
end
function Results(model::M, solution::DESolution, subsampler::SubSampler) where {M <: Model}
	SubSampledResults(model, solution, subsampler)
end

Base.iterate(r::Results, state...) = iterate(tuples(r.solution), state...)
function Base.iterate(r::SR) where {T, M <: Model{T}, SR <: SubSampledResults{M}}
	sampled = sample_space(r.solution(0), r.subsampler)
	((sampled, 0), 0)
end
function Base.iterate(r::SR, (prev_val, prev_t)) where {T, M <: Model{T}, SR <: SubSampledResults{M}}
	new_t = prev_t + r.subsampler.dt
	if new_t > r.solution.t[end]
		return nothing
	else
		sampled = sample_space(r.solution(new_t), r.subsampler)
		return ((sampled, new_t), new_t)
	end
end

function analyse(a::Analyses, results::AbstractResults{<:M}, output::AbstractOutput) where {M <: Model}
    a.plots .|> (plot_type) -> plot_and_save(plot_type, results, output)
    #return results
end

export AbstractResults, AbstractPlotSpecification, output_name
export plot_and_save
export SubSampler, sample, Analyses, Results, SubSampledResults
export spatiotemporal_data
export analyse

end
