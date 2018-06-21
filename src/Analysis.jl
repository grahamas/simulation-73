module Analysis
# * Load Modules

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
# using PyCall
# @pyimport mpl_toolkits.mplot3d as mpl
# Axes3D = mpl.Axes3D
using Plots; gr() #pyplot()

# upscale = 8 #8x upscaling in resolution
# fntsm = Plots.font("sans-serif", 10.0*upscale)
# fntlg = Plots.font("sans-serif", 14.0*upscale)
# default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
# default(size=(600*upscale,400*upscale)) #Plot canvas size
# default(dpi=300) #Only for PyPlot - presently broken

using Colors
#using PerceptualColourMaps

# * Analysis types
const Analyses = Dict

export Analyses

# * Timeseries
const PopTimeseries1D{ValueT<:Real} = Array{ValueT, 3} # 1 spatial dimension
const Timeseries1D{ValueT<:Real} = Array{ValueT,2}

# * Plotting
# ** Plot gif of solution
function plot_activity_gif(t, x, timeseries::PopTimeseries1D, safe_write_fn; file_name="solution.gif",
		      disable=0, fps=15, pop_peak_timeseries=[])
    @assert size(timeseries, 2) == 2
    if disable != 0
	return
    end
    gr()
    max_activity = maximum(timeseries, (1,2,3))[1] # I don't know why this index is here.
    min_activity = minimum(timeseries, (1,2,3))[1]
    anim = @animate for i in 1:length(t)
	Plots.plot(x, timeseries[:,1,i], lab="E",
	     ylim=(min_activity, max_activity), title="t=$(t[i])",
             xlab="Space", ylab="Proportion pop. active")
        Plots.plot!(x, timeseries[:,2,i], lab="I")
	for peak_timeseries in pop_peak_timeseries
	    Plots.scatter!(peak_timeseries[i][1], peak_timeseries[i][2], markercolor=peak_timeseries[i][3])
	end
    end
    safe_write_fn("activity.gif", (path) -> gif(anim, path, fps=floor(Int,fps)))
end

# ** Plot heatmap
using PyPlot, PyCall
@pyimport mpl_toolkits.axes_grid1 as axgrid

function plot_heatmap(t, x, timeseries, safe_write_fn; name="heatmap.png", kwargs...)
    base, ext = splitext(name)
    # @assert(all([ndims(timeseries) == 3, size(timeseries,2) == 2]))
    # heatmap(x=t, y=x, timeseries[:,1,:]; kwargs...)
    # E_file_name = "$(base)_E$ext"
    # safe_write_fn(E_file_name, savefig)
    # heatmap(x=t, y=x, timeseries[:,2,:]; kwargs...)
    # I_file_name = "$(base)_I$ext"
    # safe_write_fn(I_file_name, savefig)
    PyPlot.clf()
    ax = PyPlot.gca()
    im = PyPlot.imshow(timeseries[:,1,:], ColorMap("viridis"), origin="lower",
extent=[t[1], t[end], x[1], x[end]])
    ax[:set_ylim](x[1], x[end])
    ax[:set_xlim](t[1], t[end])
    PyPlot.xlabel("Time (s)")
    PyPlot.ylabel("Space (a.u.)")
    ax[:set_aspect](:auto)
    divider = axgrid.make_axes_locatable(ax)
    cax = divider[:append_axes]("right", size="7%", pad="2%")#, size="5%")
    PyPlot.colorbar(im, cax=cax)
    E_file_name = "$(base)_E$ext"
    safe_write_fn(E_file_name, PyPlot.savefig)

    PyPlot.clf()
    ax = PyPlot.gca()
    im = PyPlot.imshow(timeseries[:,1,:], ColorMap("viridis"), origin="lower",
                       extent=[t[1], t[end], x[1], x[end]])
    ax[:set_xlim](t[1], t[end])
    ax[:set_yticks]([])
    ax[:set_aspect](:auto)
    divider = axgrid.make_axes_locatable(ax)
    cax = divider[:append_axes]("right", size="7%", pad="2%")#, size="5%")
    PyPlot.colorbar(im, cax=cax)
    E_file_name = "$(base)_E_naked$ext"
    safe_write_fn(E_file_name, PyPlot.savefig)
end

# ** Plot solution as surface
function plot_solution_surface(solution::Timeseries1D, mesh::SpaceMesh, T, dt, safe_write_fn;
                               save=nothing, seriestype=:surface)
    time_range = 0:dt:T
    x_range = mesh.dims
    Plots.plot(time_range, x_range, solution, seriestype=seriestype)
    if save != nothing
        safe_write_fn(save, savefig)
    end
end

# ** Plot nonlinearity
function plot_nonlinearity(nonlinearity_fn, safe_write_fn, pop_names; disable=0)
    if disable != 0
        return
    end
    resolution = 100
    n_pops = length(pop_names)
    one_pop_x = linspace(-1,15,resolution)
    x_range = repeat(one_pop_x, outer=(n_pops))
    y_output = reshape(nonlinearity_fn(x_range), (:, n_pops))
    Plots.plot(one_pop_x, y_output, lab=pop_names,
         xlab="Input current", ylab="Proportion pop. reaching at least threshold")
    safe_write_fn("nonlinearity.png", savefig)
end
# * Down sampling
function down_sample(t, space::StepRangeLen, timeseries; spatial_stride=1, temporal_stride=1)
    x = x_range(mesh)
    @assert size(timeseries,1) == length(x)
    @assert size(timeseries,3) == length(t)
    spatial_stride = floor(Int, spatial_stride)
    temporal_stride = floor(Int, temporal_stride)
    t_dx = 1:temporal_stride:length(t)
    x_dx = 1:spatial_stride:length(x)
    return t[t_dx], x[x_dx], timeseries[x_dx, :, t_dx]
end
# * Run analyses

function analyse_WilsonCowan73_solution(soln, sim::Union{Simulation, ParameterSearch})
    write_fn = make_write_fn(sim.output)
    write_params(write_fn, sim)
    analyse_WilsonCowan73_solution(soln, write_fn, sim.model; sim.analyses...)
end

function analyse_WilsonCowan73_solution(soln, write_fn::Function, model::Model; down_sampling=nothing, nonlinearity=nothing,
                                        pop_names=nothing, activity_gif=nothing, heatmap=nothing)
    space = calculate(model.space)
    timeseries = soln.u
    n_pops = length(model.α)
    n_space = length(space)
    @assert all(size(timeseries, [1,2]) .== (n_space, n_pops))
    if (down_sampling != nothing) ds_t, ds_x, ds_timeseries = down_sample(soln.t, space, timeseries; down_sampling...) end
    if (heatmap != nothing) plot_heatmap(ds_t, ds_x, ds_timeseries, write_fn; heatmap...) end
    #pop_peak_timeseries = calc_pop_peak_timeseries(timeseries, 0)
    #if (nonlinearity != nothing) plot_nonlinearity(soln.prob.p.nonlinearity_fn, write_fn, pop_names; nonlinearity) end
    if (activity_gif != nothing) plot_activity_gif(ds_t, ds_x, ds_timeseries, write_fn; activity_gif...) end
end

# * Export
export analyse_WilsonCowan73_solution,

    find_local_maxima,

    plot_solution_surface

end
