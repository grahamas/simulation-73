module Analysis
using CalculatedParameters, Parameters
using Parameters
using DiffEqBase: AbstractODESolution, ODESolution, DESolution
using OrdinaryDiffEq: ODECompositeSolution, InterpolationData, CompositeInterpolationData
using Modeling
using Records

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots


# * Analysis types
abstract type AbstractPlotSpecification end
function output_name(af::AbstractPlotSpecification)
	af.output_name
end

function plot_and_save(plot_spec::APS, simulation::Simulation, output::AbstractOutput) where {APS <: AbstractPlotSpecification}
	@info "Entered plot and save"
	save_fn(fn, plt) = savefig(plt, fn)
	@info "Calling output plot"
	simulation.output(save_fn, output_name(plot_spec), plot(plot_spec, solution; plot_spec.kwargs...))
end

@with_kw struct Analyses{T}
	plots::Array{AbstractPlotSpecification}
end

function analyse(simulation::simulation)
    @info "Begin analysis."
    for plot_spec in simulation.analyses.plot_specs
    	plot_and_save(plot_spec, simulation)
    end
end

export AbstractPlotSpecification, output_name
export plot_and_save
export Analyses, analyse

end
