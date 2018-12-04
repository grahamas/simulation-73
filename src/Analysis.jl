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
abstract type AbstractFigure end
function output_name end

function plot_and_save(plot_obj::AF, results::AbstractResults, output::Output) where {AF <: AbstractFigure}
	save_fn(fn, plt) = savefig(plt, fn)
	output(save_fn, "$(output_name(plot_obj)).png", plot(plot_obj, results; plot_obj.kwargs...))
end


@with_kw struct SubSampler
	dt::Float64
	space_strides::Array{Integer}
end

function sample end

@with_kw struct Analyses
	subsampler::Union{SubSampler,Nothing} = nothing
	plots::Array{AbstractFigure}
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

function Results(model::M, solution::DESolution, subsampler::Nothing) where {M <: Model}
	Results(model, solution)
end
function Results(model::M, solution::DESolution, subsampler::SubSampler) where {M <: Model}
	SubSampledResults(model, solution, subsampler)
end

function spatiotemporal_data(r::Results)
	return spatiotemporal_data(r.solution, r.model)
end

function spatiotemporal_data(r::SubSampledResults)
	return sample(r.subsampler, r.solution, r.model)
end

function analyse(a::Analyses, results::AbstractResults{<:M}, output::Output) where {M <: Model}
    a.plots .|> (plot_type) -> plot_and_save(plot_type, results, output)
    #return results
end

export AbstractResults, AbstractFigure, output_name
export plot_and_save
export SubSampler, sample, Analyses, Results, SubSampledResults
export spatiotemporal_data
export analyse

end
