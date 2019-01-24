module Analysis
using CalculatedParameters, Parameters
using Parameters
using DiffEqBase: AbstractODESolution, ODESolution, DESolution
using OrdinaryDiffEq: ODECompositeSolution, InterpolationData, CompositeInterpolationData
using Modeling
using Records
using Simulating

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots

# * Analysis types are defined in Simulating.jl

function output_name(af::AbstractPlotSpecification)
	af.output_name
end

function plot_and_save(plot_spec::APS, simulation::Simulation) where {APS <: AbstractPlotSpecification}
	@info "Entered plot and save"
	save_fn(fn, plt) = savefig(plt, fn)
	@info "Calling output plot"
	plot_obj = plot(plot_spec, simulation; plot_spec.kwargs...)
	simulation.output(save_fn, output_name(plot_spec), plot_obj)
end

function analyse(simulation::Simulation)
    @info "Begin analysis."
    for plot_spec in simulation.analyses.plot_specs
    	plot_and_save(plot_spec, simulation)
    end
end

export output_name
export plot_and_save
export analyse

end
