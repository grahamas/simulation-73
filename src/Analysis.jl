module Analysis
using CalculatedParameters
using Parameters
using DiffEqBase
using Modeling
using Records


# * Analysis types

abstract type AbstractResults{M<:Model} end
abstract type AbstractFigure end
struct PlotResults{AF <: AbstractFigure}
    results::AbstractResults
end

function plot_and_save(plot_obj::AF, base_name::AbstractString, results::AbstractResults, output::Output) where {AF <: AbstractFigure}
	save_fn(fn, plt) = save(plt, fn)
	output(save_fn, base_name, plot(PlotResults{AF}(results); plot_obj.kwargs...))
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

function (a::Analyses{M})(results::Results{M}, output::Output) where {M <: Model}
    a.plots .|> (plot_fn) -> plot_fn(results)
end

export analyse, Analyses, AbstractFigure, SubSampler, spatiotemporal_data,
		AbstractResults, Results, SubSampledResults, plot_and_save,
		PlotResults, sample

end
