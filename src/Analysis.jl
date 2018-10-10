module Analysis
using CalculatedParameters
using Parameters
using DiffEqBase

# * Analysis types
@with_kw struct Analyses{M <: Model}
	subsampler::V where {V <: Union{SubSampler{M},Nothing}}
	plots::Array{AbstractPlot}
end

(a::Analyses{M})(results::Results{M}) where {M <: Model}
    a.plots .|> (plot_fn) -> plot_fn(results)
}
end

abstract type AbstractResults{M<:Model} end
@with_kw struct Results{M} <: AbstractResults{M} where {M <: Model}
	model::M
	soln::DESolution
end
@with_kw struct SubSampledResults{M} <: AbstractResults{M} where {M <: Model}
	model::M
	soln::DESolution
	subsampler::SubSampler{M}
end

function spatiotemporal_data(r::Results)
	return spatiotemporal_data(r.soln, r.model)
end

function spatiotemporal_data(r::SubSampledResults)
	return r.subsampler(r.soln, r.model)
end

function spatiotemporal_data(soln::DESolution, model::Model)
	t = soln.t
	x = space(model)
	u = soln.u
	return (t,x,u)
end

struct SubSampler{M <: Model}
	dt::Float64
	space_strides::Array{Integer}
end

abstract type AbstractPlot end

export analyse, Analyses

end
