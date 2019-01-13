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
using WCMNonlinearity

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
    x = space_arr(simulation)
    t = time_arr(simulation)
    max_val = maximum(simulation)
    @animate for time_dx in 1:length(t) # TODO @views
        plot(x, solution[:, 1, time_dx]; label=pop_names[1],
            ylim=(0,max_val), title="t = $(round(t[time_dx], digits=4))", kwargs...)
        for i_pop in 2:length(pop_names)
            plot!(x, solution[:, i_pop, time_dx]; label=pop_names[i_pop], kwargs...)
        end

    end
end

struct SpaceTimePlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
SpaceTimePlot(; output_name = "spacetime.png", kwargs...) = SpaceTimePlot(output_name, kwargs)
@recipe function f(plot_spec::SpaceTimePlot, simulation::Simulation{T,M}) where {T,M<:WCMSpatial1D}
    v_space = space_arr(simulation)
    v_time = time_arr(simulation)
    clims := (minimum(simulation), maximum(simulation))
    grid := false
    layout := (2,1)
    for i_pop in 1:2 # TODO!!
        @series begin
            seriestype --> :heatmap
            subplot := i_pop
            x := v_time
            y := v_space
            simulation.solution[:,i_pop,:]
        end
    end
end

export SpaceTimePlot

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

    one_pop_x = 
    #delete!.(Ref(plotattributes),[:resolution,:fn_bounds])

    xlab := "Input current"
    ylab := "Proportion pop. reaching at least threshold"

    one_pop_x = range(fn_bounds[1], stop=fn_bounds[2], length=resolution)

    for i_pop in 1:length(pop_names)
        @series begin
            lab --> pop_names[i_pop]
            seriestype := :line
            x := one_pop_x
            y := nonlinearity(nonlinearity_fns[i_pop], collect(one_pop_x))
            ()
        end
    end
end

export NonlinearityPlot


struct NeumanTravelingWavePlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
NeumanTravelingWavePlot(; output_name="traveling_wave.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = NeumanTravelingWavePlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::NeumanTravelingWavePlot{T}, simulation::Simulation{T,M}) where {T,M<:WCMSpatial1D}
    t = time_arr(simulation)
    space_origin::Int = get_origin(simulation) # TODO: Remove 1D return assumption
    di = max(1, round(Int, simulation.solver.simulated_dt / plot_spec.dt))
    x = space_arr(simulation)[space_origin:di:end] # TODO: Remove 1D return assumption
    for time_dx in 1:length(t)
        @series begin
            seriestype := :line
            x := x
            y := simulation.solution[space_origin:di:end,:,time_dx] * [1.0, -1.0] # Subtract inhibitory activity...
            ()
        end
    end
end

export NeumanTravelingWavePlot #, TravelingWavePlot

end