module WCMAnalysis

using Modeling
using Analysis
using Records
import DifferentialEquations: DESolution
using Meshes
using Colors
using RecipesBase
using Parameters

@with_kw struct Animate <: AbstractFigure
    kwargs::Dict = Dict()
    fps::Int = 20
end
Animate(; fps, kwargs...) = Animate(; fps=fps, kwargs = kwargs)
function Analysis.plot_and_save(plot_type::Animate, base_name::AbstractString, results::AbstractResults, output::Output)
    save_fn(fn, anim) = mp4(plt, fn; fps=plot_type.fps)
    output(save_fn, base_name, plot(PlotResults{Animate}(results); plot_type.kwargs...))
end
# TODO: Implement animation (using RecipesBase being the challenge...)
# NOTE: Probably requires PR to Plots.jl

@with_kw struct SpaceTimePlot <: AbstractFigure
    kwargs::Dict = Dict()
end
SpaceTimePlot(; kwargs...) = SpaceTimePlot(; kwargs = kwargs)
@recipe function f(h::PlotResults{SpaceTimePlot})
    v_time, v_space, timeseries = spatiotemporal_data(h.results)
    @assert size(timeseries, 2) == 2      # only defined for 2 pops
    clims := (minimum(timeseries), maximum(timeseries))
    grid := false
    layout := (2,1)
    for i_pop in 1:size(timeseries,2)
        @series begin
            seriestype --> :heatmap
            subplot := i_pop
            x := v_time
            y := v_space
            timeseries[:,i_pop,:]
        end
    end
end

export SpaceTimePlot

# ** Plot nonlinearity
@with_kw struct NonlinearityPlot <: AbstractFigure
    kwargs::Dict = Dict()
end
NonlinearityPlot(; kwargs...) = NonlinearityPlot(; kwargs = kwargs)
@recipe function f(n::PlotResults{NonlinearityPlot}; resolution=100, fn_bounds=(-1,15))
    pop_names = n.results.model.pop_names
    nonlinearity_fns = get_value.(Calculated(n.results.model).nonlinearity)
    n_pops = length(pop_names)

    one_pop_x = linspace(fn_bounds..., resolution)
    delete!.(plotattributes,[:resolution,:fn_bounds])

    xlab := "Input current"
    ylab := "Proportion pop. reaching at least threshold"

    for i_pop in 1:length(pop_names)
        @series begin
            lab --> pop_names[i_pop]
            seriestype := :line
            x := one_pop_x
            y := nonlinearity_fns[i_pop](one_pop_x)
            ()
        end
    end
end

export NonlinearityPlot

end