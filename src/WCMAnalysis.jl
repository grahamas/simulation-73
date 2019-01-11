module WCMAnalysis

using Modeling
using Analysis
using Records
import DifferentialEquations: DESolution
using Meshes
using RecipesBase
using Plots
using Parameters
using Memoize
using WCM
using CalculatedParameters
using Simulating

struct Animate <: AbstractPlotSpecification
    fps::Int
    output_name::String
    kwargs::Dict
end
Animate(; fps=20, output_name="animation.mp4", kwargs...) = Animate(fps, output_name, kwargs)
function Analysis.plot_and_save(plot_spec::Animate, simulation::Simulation)
    save_fn(name, anim) = mp4(anim, name; fps=plot_spec.fps)
    simulation.output(save_fn, output_name(plot_spec), animate(simulation; plot_spec.kwargs...))
end

export Animate

function RecipesBase.animate(simulation::Simulation{T,M}; kwargs...) where {T,M<:WCMSpatial1D}
    solution = simulation.solution
    pop_names = simulation.model.pop_names
    x = get_space_arr(simulation)
    t = get_time_arr(simulation)
    @animate for time_dx in 1:length(t) # TODO @views
        plot(x, solution[:, 1, time_dx]; label=pop_names[1],
            ylim=(0,max_val), title="t = $(round(t[time_dx], digits=4))", kwargs...)
        for i_pop in 2:length(pop_names)
            plot!(x, solution[:, i_pop, time_dx]; label=pop_names[i_pop], kwargs...)
        end

    end
end


# struct SpaceTimePlot <: AbstractPlotSpecification
#     output_name::String
#     kwargs::Dict
# end
# SpaceTimePlot(; output_name = "spacetime.png", kwargs...) = SpaceTimePlot(output_name, kwargs)
# @recipe function f(plot_spec::SpaceTimePlot, results::AbstractResults; kwargs...)
#     v_space = get_space(results)
#     v_time = get_time(results)
#     clims := (unsampled_minimum(results), unsampled_maximum(results))
#     grid := false
#     layout := (2,1)
#     for i_pop in 1:size(timeseries,2)
#         @series begin
#             seriestype --> :heatmap
#             subplot := i_pop
#             x := v_time
#             y := v_space
#             get_pop(results, i_pop)
#         end
#     end
# end

# export SpaceTimePlot

# ** Plot nonlinearity
struct NonlinearityPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
NonlinearityPlot(; output_name = "nonlinearity.png", kwargs...) = NonlinearityPlot(output_name, kwargs)

@recipe function f(plot_spec::NonlinearityPlot, simulation::Simulation{T,M}; resolution=100, fn_bounds=(-1.0,15.0)) where {T,M<:WCMSpatial1D}
    pop_names = simulation.model.pop_names
    n_pops = length(pop_names)

    nonlinearity_fns = get_value.(Calculated(simulation.model).nonlinearity)

    one_pop_x = range(fn_bounds[1], stop=fn_bounds[2], length=resolution)
    #delete!.(Ref(plotattributes),[:resolution,:fn_bounds])

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

# struct TravelingWavePlot{T} <: AbstractPlotSpecification
#     output_name::String
#     dt::T
#     kwargs::Dict
# end
# TravelingWavePlot{T}(; output_name="traveling_wave.png", dt=nothing, kwargs...) where T = TravelingWavePlot{T}(output_name, dt, kwargs)
# @recipe function f(plot_spec::TravelingWavePlot, results::AbstractResults{WCMSpatial1D{T,C,N,S}}; kwargs...) where {T,C,N,S}
#     sampled_results = resample(results, dt=plot_spec.dt)
#     space = get_space(results)
#     for (frame, t) in sampled_results
#         @series begin
#             seriestype := :line
#             x := space
#             y := frame
#             ()  # TOOOOOOOOOOODOOOOOOOOOOOOOOOOOOOOOOOO
#         end
#     end
# end

struct NeumanTravelingWavePlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
NeumanTravelingWavePlot(; output_name="traveling_wave.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = NeumanTravelingWavePlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::NeumanTravelingWavePlot{T}, simulation::Simulation{T,M}) where {T,M<:WCMSpatial1D}
    @info "entered plot"
    t = get_time_arr(simulation)
    space = get_space_arr(simulation)
    space_origin = findfirst((x) -> x ≈ 0.0, space)
    @info "looping"
    for time_dx in 1:length(t)
        @info "loop $(t[time_dx])"
        @series begin
            seriestype := :line
            x := spaceF
            y := frame[space_origin:end,:,time_dx] * [1.0, -1.0] # Subtract inhibitory activity...
            ()
        end
    end
end

export NeumanTravelingWavePlot #, TravelingWavePlot

end