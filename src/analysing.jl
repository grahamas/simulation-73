abstract type AbstractAnalysis end

abstract type AbstractPlotSpecification <: AbstractAnalysis end
abstract type AbstractSpaceTimePlotSpecification <: AbstractPlotSpecification end

# Default that should work for most
function output_name(aps::AbstractPlotSpecification)
	aps.output_name
end

function plot_and_save(plot_spec::AbstractPlotSpecification, execution::Execution, output_dir)
	plot_obj = plot(plot_spec, execution; plot_spec.kwargs...)
	safesave(joinpath(output_dir, output_name(plot_spec)), plot_obj)
end

function analyse(plot_spec::AbstractPlotSpecification, execution::Execution, output_dir)
   	plot_and_save(plot_spec, execution, output_dir)
end
