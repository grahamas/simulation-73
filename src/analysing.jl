abstract type AbstractAnalysis end

abstract type AbstractPlotSpecification <: AbstractAnalysis end
abstract type AbstractSpaceTimePlotSpecification <: AbstractPlotSpecification end

# Default that should work for most
function output_name(aps::AbstractPlotSpecification)
	aps.output_name
end

function plot_and_save(plot_spec::AbstractPlotSpecification, execution::Execution, output_dir::AbstractString)
	path = joinpath(output_dir, output_name(plot_spec))
	DrWatson.recursively_clear_path(path)
	plot_obj = plot(plot_spec, execution; plot_spec.kwargs...)
	savefig(plot_obj, path)
end

function analyse(plot_spec::AbstractPlotSpecification, execution::Execution, output_dir)
   	plot_and_save(plot_spec, execution, output_dir)
end
