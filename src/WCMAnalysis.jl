module WCMAnalysis

using Modeling
using Analysis
using Records
import DifferentialEquations: DESolution
using Meshes
using Colors
using RecipesBase
using Memoize

(a::Analyses{WCMSpatial1D})(results::Results{WCMSpatial1D})
    a.plots .|> (plot_fn) -> plot_fn(results)
end

# * Plotting
macro defined_userplot(typename::Symbol)
    @shorthands((Symbol ∘ lowercase ∘ string)(typename))
end

const PopActivity1D = Array{Float64, 3}
@userplot struct SpaceTimePlot
    results::AbstractResults
end
@recipe function f(h::SpaceTimePlot)
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

# ** Plot nonlinearity
@userplot struct NonlinearityPlot
    results::AbstractResults
end
@recipe function f(n::NonlinearityPlot; resolution=100, fn_bounds=(-1,15))
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

function SubSampler{WCMSpatial1D}(; dt::Float64=0, spatial_stride::Int=1)
    @assert dt > 0
    SubSampler{WCMSpatial1D}(dt, [spatial_stride])
end

@memoize
(subsampler::SubSampler{M})(soln::DESolution, model::M) where {M <: WCMSpatial1D}
    timepoints = minimum(soln.t):subsampler.dt:maximum(soln.t)
    # Assuming densely sampled.
    timesampled_u = soln(timepoints)

    space = space(model)
    space_stride = subsampler.space_strides[1]
    sampled_space = space[1:space_stride:end]

    sampled_u = timesampled_u[1:space_stride:end,:,:]

    return timepoints, sampled_space, sampled_u
end

end