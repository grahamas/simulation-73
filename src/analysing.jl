abstract type AbstractAnalysis end

abstract type AbstractPlotSpecification <: AbstractAnalysis end
abstract type AbstractSpaceTimePlotSpecification <: AbstractPlotSpecification end

# Default that should work for most
function output_name(aps::AbstractPlotSpecification)
	aps.output_name
end

function plot_and_save(plot_spec::AbstractPlotSpecification, simulation::Simulation, output::AbstractOutput)
	save_fn(fn, plt) = savefig(plt, fn)
	plot_obj = plot(plot_spec, simulation; plot_spec.kwargs...)
	output(save_fn, output_name(plot_spec), plot_obj)
end

function analyse(plot_spec::AbstractPlotSpecification, simulation::Simulation, output::AbstractOutput)
   	plot_and_save(plot_spec, simulation, output)
end