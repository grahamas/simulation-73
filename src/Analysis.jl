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


@with_kw struct SubSampler{M <: Model}
	dt::Float64
	space_strides::Array{Integer}
end

function sample end

@with_kw struct Analyses{M <: Model}
	subsampler::V where {V <: Union{SubSampler{M},Nothing}}
	plots::Array{AbstractFigure}
end

@with_kw struct Results{M} <: AbstractResults{M}
	model::M
	soln::DESolution
end
@with_kw struct SubSampledResults{M} <: AbstractResults{M}
	model::M
	soln::DESolution
	subsampler::SubSampler{M}
end

function spatiotemporal_data(r::Results)
	return spatiotemporal_data(r.soln, r.model)
end

function spatiotemporal_data(r::SubSampledResults)
	return sample(r.subsampler, r.soln, r.model)
end

function spatiotemporal_data(soln::DESolution, model::Model)
	t = soln.t
	x = space(model)
	u = soln.u
	return (t,x,u)
end

function analyse(a::Analyses{M}, results::Results{<:M}, output::Output) where {M <: Model}
    a.plots .|> (plot_type) -> plot_and_save(plot_type, results, output)
end

export AbstractResults, AbstractFigure, output_name
export plot_and_save
export SubSampler, sample, Analyses, Results, SubSampledResults
export spatiotemporal_data
export analyse

end
