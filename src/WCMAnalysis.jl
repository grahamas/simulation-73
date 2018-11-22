module WCMAnalysis

using Modeling
using Analysis
using Records
import DifferentialEquations: DESolution
using Meshes
using Colors
using RecipesBase
using Plots; pyplot()
using Parameters
using Memoize
using WCM
using CalculatedParameters

struct Animate <: AbstractFigure
    fps::Int
    kwargs::Dict
end
Animate(; fps=20, kwargs...) = Animate(fps, kwargs)
function Analysis.plot_and_save(plot_type::Animate, results::AbstractResults, output::Output)
    save_fn(name, anim) = mp4(anim, name; fps=plot_type.fps)
    output(save_fn, "animation.mp4", animate(results; plot_type.kwargs...))
end
# TODO: Implement animation (using RecipesBase being the challenge...)
# NOTE: Probably requires PR to Plots.jl

export Animate

struct SpaceTimePlot <: AbstractFigure
    kwargs::Dict
end
SpaceTimePlot(; kwargs...) = SpaceTimePlot(kwargs)
@recipe function f(plot_type::SpaceTimePlot, results::AbstractResults)
    v_time, v_space, timeseries = spatiotemporal_data(results)
    timeseries = cat(timeseries..., dims=3)
    @assert (size(timeseries, 2) == 2) size(timeseries)   # only defined for 2 pops
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
Analysis.output_name(plt::SpaceTimePlot) = "spacetimeplot"

export SpaceTimePlot

# ** Plot nonlinearity
struct NonlinearityPlot <: AbstractFigure
    kwargs::Dict
end
NonlinearityPlot(; kwargs...) = NonlinearityPlot(kwargs)
@recipe function f(plot_type::NonlinearityPlot, results::AbstractResults; resolution=100, fn_bounds=(-1,15))
    pop_names = results.model.pop_names
    nonlinearity_fns = get_value.(Calculated(results.model).nonlinearity)
    n_pops = length(pop_names)

    one_pop_x = range(fn_bounds[1], stop=fn_bounds[2], length=resolution)
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
Analysis.output_name(plt::NonlinearityPlot) = "nonlinearityplot"

export NonlinearityPlot

function Analysis.SubSampler(dt::Float64, spatial_stride::Int) where WCMSpatial1D
    @assert dt > 0
    SubSampler{WCMSpatial1D}(dt, [spatial_stride])
end

@memoize function Analysis.sample(subsampler::SubSampler{WCMSpatial1D}, soln::DESolution, model::WCMSpatial1D)
    println("Sampling WCMSpatial1D")
    timepoints = minimum(soln.t):subsampler.dt:maximum(soln.t)
    # Assuming densely sampled.
    timesampled_u = soln(timepoints)

    space = space_array(model)
    space_stride = subsampler.space_strides[1]
    sampled_space = space[1:space_stride:end]

    sampled_u = timesampled_u[1:space_stride:end,:,:]
    @show size(sampled_u)

    return timepoints, sampled_space, sampled_u
end

function Analysis.spatiotemporal_data(soln::DESolution, model::WCMSpatial1D)
    t = soln.t
    x = space_array(model)
    u = soln.u
    return (t,x,u)
end

function RecipesBase.animate(results::AbstractResults{WCMSpatial1D}; kwargs...)
    t, x, data = spatiotemporal_data(results)
    pop_names = results.model.pop_names
    @assert length(x) == size(data,1)
    @animate for i_time in 1:length(t)
        plot(x, data[:, 1, i_time]; label=pop_names[1], kwargs...)
        for i_pop in 2:size(data,2)
            plot!(x, data[:, 1, i_time]; label=pop_names[i_pop], kwargs...)
        end

    end
    # Not using animate(results.solution) to use subsampling
end

end