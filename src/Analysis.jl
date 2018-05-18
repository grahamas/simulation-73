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
using PerceptualColourMaps

# * Timeseries
const PopTimeseries1D{ValueT<:Real} = Array{ValueT, 3} # 1 spatial dimension
const Timeseries1D{ValueT<:Real} = Array{ValueT,2}

# * Moving peaks
# ** Definition
struct MovingPeak{ValueT<:Real}
    heights::Array{ValueT}
    position_indices::Array{Int}
    time_indices::Array{Int}
    scions::Array{MovingPeak{ValueT}}# = MovingPeak[]
    sources::Array{MovingPeak{ValueT}}# = MovingPeak[] #Can be more than one if collision
end

# ** Constructor
MovingPeak(heights, locations, times) = MovingPeak(heights, locations, times,
						   MovingPeak[], MovingPeak[])

# ** Methods
cur_time(peak::MovingPeak) = peak.time_indices[end]
add_scions!(peak::MovingPeak, new_scions) = append!(peak.scions, new_scions)
function add_source!(peak::MovingPeak, new_source)
    @assert peak != new_source "Tried to make peak its own source."
    if new_source ∉ peak.sources
	push!(peak.sources, new_source)
    end
end
function pos_vel_acc(peak::MovingPeak)
    return (pos(peak), vel(peak), acc(peak))
end

function pos(peak::MovingPeak, offset::Int=0)
    if (length(peak.position_indices) + offset) < 1
	if length(peak.sources) > 0
	    return maximum(pos.(peak.sources, offset+1))
	else
	    return NaN
	end
    else
	return peak.position_indices[end+offset]
    end
end
vel(peak::MovingPeak, offset::Int=0) = pos(peak,offset) - pos(peak,offset-1)
acc(peak::MovingPeak) = vel(peak) - vel(peak,-1)

  function pos_vel_acc(peak::MovingPeak, putative_continuation::MovingPeak)
      put_pos = pos(putative_continuation)
      put_vel = put_pos - pos(peak)
      put_acc = put_vel - vel(peak)
      return (put_pos, put_vel, put_acc)
  end

# * Find moving peaks
function find_moving_peaks(timeseries::Timeseries1D{ValueT},rtol)::Array{MovingPeak} where {ValueT <: Real}
    peaks = MovingPeak[]
    n_time = size(timeseries,2)
    for i_time = 1:n_time
	space_frame = timeseries[:,i_time]
	new_maxima = find_local_maxima(space_frame)
	new_peaks = [MovingPeak(ValueT[space_frame[idx]], Int[idx], Int[i_time])
		     for idx in new_maxima]
	continue_moving_peaks!(peaks, new_peaks, rtol)
    end
    return peaks
end

function find_moving_peaks(timeseries::PopTimeseries1D, rtol)
    [find_moving_peaks(timeseries[:,j,:], rtol) for j in 1:size(timeseries,2)]
end
# from: https://discourse.julialang.org/t/how-to-identify-local-maxima-peaks-in-a-time-signal/6000/2
function find_local_maxima(signal::Vector, threshold=1e-3)
    inds = Int[]
    if length(signal)>1
	if signal[1]>signal[2] && signal[1] > threshold
	    push!(inds,1)
	end
	for i=2:(length(signal)-1)
	    if signal[i-1]<signal[i]>signal[i+1] && signal[i] > threshold
		push!(inds,i)
	    end
	end
	if signal[end]>signal[end-1] && signal[end] > threshold
	    push!(inds,length(signal))
	end
    end
    inds
end
function accept_all_potential_continuations!(peaks::Array{MovingPeak},
					   next_peaks::Array{MovingPeak},
					   continuations::Array{Array{Int,1},1},
					   sources::Array{Array{Int,1},1})
    eaten_next_peaks = Set{Int}()
    for (i_pre, pre_peak) in enumerate(peaks)
	if length(continuations[i_pre]) == 1
	    continuation_dx = continuations[i_pre][1]
	    if length(sources[continuation_dx]) == 1
		continue_peak!(pre_peak, next_peaks[continuation_dx])
		push!(eaten_next_peaks, continuation_dx)
		continue
	    end
	end
    end
    for (i_pre, pre_peak) in enumerate(peaks)
	# No one-to-one mapping
	scions = [next_peaks[i] for i ∈ continuations[i_pre] if i ∉ eaten_next_peaks]
	add_scions!(pre_peak, scions)
	for scion in scions
	    add_source!(scion, pre_peak)
	end
    end
    deleteat!(next_peaks, sort!(collect(eaten_next_peaks)))
    append!(peaks, next_peaks)
end
function find_potential_continuations(peaks::Array{MovingPeak},
				      next_peaks::Array{MovingPeak}, rtol)
    # Array of array of indices pointing to potential continuations in
    # next_peaks. Probably needs to be initialized.
    potential_continuations = Array{Int,1}[[] for _ in peaks]
    # Same, but sources in peaks.
    potential_sources = Array{Int,1}[[] for _ in next_peaks]
    for (i_peak, peak) in enumerate(peaks)
	for (i_next_peak, next_peak) in enumerate(next_peaks)
	    if is_continuation(peak, next_peak, rtol)
		push!(potential_continuations[i_peak], i_next_peak)
		push!(potential_sources[i_next_peak], i_peak)
	    end
	end
    end
    return (potential_continuations, potential_sources)
end

function is_continuation(peak::MovingPeak, next_peak::MovingPeak,
			 rtol=0, atol=1) #permissive arbitrary??
    # should depend on dt
    if cur_time(peak) != cur_time(next_peak) - 1
	return false
    end
    cur_pos, cur_vel, cur_acc = pos_vel_acc(peak)
    if isnan(cur_acc)
	cur_acc = 0
    end
    if isnan(cur_vel)
        cur_vel = 0
    end
    proj_vel = cur_vel + cur_acc
    proj_pos = cur_pos + proj_vel
    next_pos, next_vel, next_acc = pos_vel_acc(peak, next_peak)
    return ≈((proj_pos - next_pos), proj_vel, rtol=rtol, atol=atol)
end
function refine_continuations!(peaks::Array{MovingPeak})
  return
end

function continue_moving_peaks!(peaks::Array{MovingPeak},
				next_peaks::Array{MovingPeak}, rtol)
    # figure out if the "next_peaks" are new peaks or continuations of existing peaks
    potential_continuations, potential_sources = find_potential_continuations(peaks, next_peaks, rtol)
    accept_all_potential_continuations!(peaks, next_peaks, potential_continuations, potential_sources)
    refine_continuations!(peaks)
end
function continue_peak!(peak::MovingPeak, next_peak::MovingPeak)
    @assert next_peak.time_indices[1] > peak.time_indices[end]
    @assert length(next_peak.sources) == 0
    append!(peak.time_indices, next_peak.time_indices)
    append!(peak.position_indices, next_peak.position_indices)
    append!(peak.heights, next_peak.heights)
end

# * Convert peaks to timeseries
function calc_pop_peak_timeseries(timeseries::PopTimeseries1D, rtol)
    pop_peaks =  find_moving_peaks(timeseries, rtol)
    n_colors = maximum(length.(pop_peaks))
    end_t_dx = size(timeseries,3)
    return [peaks_to_timeseries(pop_peaks[j], end_t_dx, n_colors) for j in 1:size(timeseries,2)]
end

function peaks_to_timeseries(peaks::Array{MovingPeak{ValueT}}, end_time_dx::Int, n_colors::Int=0) where {ValueT <: Real}
    if n_colors > 0
	colors = distinguishable_colors(n_colors)
	timeseries = Tuple{Array{Int}, Array{ValueT}, Array{RGBA}}[(Int[], ValueT[], RGBA[])
						       for _ in 1:end_time_dx]
    else
	timeseries = Tuple{Array{Int}, Array{ValueT}}[(Int[], ValueT[]) for _ in 1:end_time_dx]
    end
    for (i_peak, peak) in enumerate(peaks)
	for (i_time_indices, time_dx) in enumerate(peak.time_indices)
	    push!(timeseries[time_dx][1], peak.position_indices[i_time_indices])
	    push!(timeseries[time_dx][2], peak.heights[i_time_indices])
	    push!(timeseries[time_dx][3], colors[i_peak])
	end
    end
    return timeseries
end

# * File ops
function make_individual_output_folder(root, simulation_name, mod_name)
    now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
    if length(mod_name) > 0
        now = join([mod_name, now], "_")
    end
    dir_name = joinpath(root, simulation_name, now)
    mkpath(dir_name)
    return dir_name
end

function make_experiment_output_folder(root, simulation_name, mod_name, experiment_name)
    @assert length(mod_name) > 0
    #now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
    #experiment_dir_name = join([experiment_name, now], "_")
    dir_name = joinpath(root, experiment_name)
    mkpath(dir_name)
    return (dir_name, if (length(simulation_name) > 0) join([simulation_name, mod_name], "_") else mod_name end)
end

function make_write_fn(; root="", simulation_name="", mod_name="", experiment_name="", other...)
    @assert(all(length.([simulation_name, root]) .> 0))
    if length(experiment_name) == 0
        dir_name, prefix = (make_individual_output_folder(root, simulation_name, mod_name), "")
    else
        dir_name, prefix = (make_experiment_output_folder(root, simulation_name, mod_name, experiment_name))
    end
    function safe_write_fn(base_name, write_fn)
        prefixed_name = join([prefix, base_name], "_")
        full_path = joinpath(dir_name, prefixed_name)
        if !(isfile(full_path))
   	     write_fn(full_path)
	else
	     warn("Tried to write existing file: $full_path")
	end
    end
    return safe_write_fn
end

function write_params(safe_write_fn; params...)
    base_name = "parameters.json"
    safe_write_fn(base_name, (path) -> write(path, JSON.json(params,4)))
end

# * Unflatten timeseries
function standardize_timeseries(timeseries, mesh::M)::PopTimeseries1D where M <: AbstractMesh
    # Join array of arrays into matrix Other Dims x Time
    cat(true_ndims(mesh)+1, [standardize_frame(frame, mesh) for frame in timeseries]...)
end
function standardize_frame(frame, mesh::FlatMesh)
    reshape(frame, size(mesh.pop_mesh))
end
function standardize_frame(frame, mesh::PopMesh)
    frame # The PopMesh shape is the standard.
end

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
function down_sample(t, mesh, timeseries; spatial_stride=1, temporal_stride=1)
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

function analyse_WilsonCowan73_solution(soln; output=nothing, analyses=nothing, other_params...)
    @assert analyses != nothing
    write_fn = make_write_fn(; output...)
    write_params(write_fn; output=output, analyses=analyses, other_params...)
    analyse_WilsonCowan73_solution(soln, write_fn; analyses...)
end

function analyse_WilsonCowan73_solution(soln, write_fn::Function; down_sampling=nothing, nonlinearity=nothing,
                                        pop_names=nothing, activity_gif=nothing, heatmap=nothing)
    timeseries = standardize_timeseries(soln.u, soln.prob.p.mesh)
    if (down_sampling != nothing) ds_t, ds_x, ds_timeseries = down_sample(soln.t, soln.prob.p.mesh, timeseries; down_sampling...) end
    if (heatmap != nothing) plot_heatmap(ds_t, ds_x, ds_timeseries, write_fn; heatmap...) end
    #pop_peak_timeseries = calc_pop_peak_timeseries(timeseries, 0)
    #if (nonlinearity != nothing) plot_nonlinearity(soln.prob.p.nonlinearity_fn, write_fn, pop_names; nonlinearity) end
    if (activity_gif != nothing) plot_activity_gif(ds_t, ds_x, ds_timeseries, write_fn; activity_gif...) end
end

# * Export
export analyse_WilsonCowan73_solution,

    find_local_maxima,

    plot_solution_surface
