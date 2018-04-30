# * Load Modules

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using PyCall
@pyimport mpl_toolkits.mplot3d as mpl
Axes3D = mpl.Axes3D
using Plots; pyplot()

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
macro safe_write(path, writer)
    quote
	if !(isfile($(esc(path))))
	    $(@schedule esc(writer))
	else
	    warn("Tried to write existing file: $(esc(path))")
	end
    end
end

function output_dir_name(; root=nothing, simulation_name=nothing, other...)
    now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
    dir_name = joinpath(root, simulation_name, now)
    mkpath(dir_name)
    return dir_name
end

function write_params(dir_name; params...)
    save_path = joinpath(dir_name, "parameters.json")
    @safe_write(save_path, write(save_path, JSON.json(params,4)))
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
function solution_gif(t, timeseries::PopTimeseries1D; dir_name="", file_name="solution.gif",
		      disable=0, subsample=1, fps=15, pop_peak_timeseries=[],
		      spatial_subsample_to=0)
    @assert size(timeseries, 2) == 2
    if disable != 0
	return
    end

    if spatial_subsample_to > 0
	spatial_stride = round(Int, size(timeseries,1) / spatial_subsample_to)
    else
	spatial_stride = 1
    end
    max_activity = maximum(timeseries, (1,2,3))[1] # I don't know why this index is here.
    min_activity = minimum(timeseries, (1,2,3))[1]
    subsample = floor(Int, subsample)
    indices = 1:spatial_stride:size(timeseries,1)
    anim = @animate for i in 1:subsample:length(t)
	plot([indices, indices], [timeseries[1:spatial_stride:end,1,i], timeseries[1:spatial_stride:end,2,i]],
	     ylim=(min_activity, max_activity), title="t=$(t[i])", legend=:none)
	for peak_timeseries in pop_peak_timeseries
	    scatter!(peak_timeseries[i][1], peak_timeseries[i][2], markercolor=peak_timeseries[i][3])
	end
    end
    save_path = joinpath(dir_name, file_name)
    @safe_write(save_path, gif(anim, save_path, fps=floor(Int,fps)))
end
# ** Plot solution as surface
function plot_solution_surface(solution::Timeseries1D, mesh::SpaceMesh, T, dt; save=nothing, seriestype=:surface)
    time_range = 0:dt:T
    x_range = mesh.dims
    plot(time_range, x_range, solution, seriestype=seriestype)
    if save != nothing
        @safe_write(save, savefig(save))
    end
end

# ** Plot nonlinearity
function plot_nonlinearity(nonlinearity_fn; dir_name=nothing, model=nothing, pop_names=nothing, other_params...)
    @assert all([dir_name, model, pop_names] .!= nothing)
    x_range = -1:0.01:15
    y_output = nonlinearity_fn(x_range)
    plot(x_range, y_output, lab=Array{String,1}(pop_names))
    savefig(joinpath(dir_name,"nonlinearity.png"))
end
# * Run analyses
function analyse_WilsonCowan73_solution(soln; analyses=nothing, other_params...)
    dir_name = output_dir_name(; analyses...)
    write_params(dir_name; analyses=analyses, other_params...)
    timeseries = standardize_timeseries(soln.u, soln.prob.p.mesh)
    #pop_peak_timeseries = calc_pop_peak_timeseries(timeseries, 0)
    plot_nonlinearity(soln.prob.p.nonlinearity_fn; dir_name=dir_name, analyses..., other_params...)
    solution_gif(soln.t, timeseries; dir_name=dir_name, #pop_peak_timeseries=pop_peak_timeseries,
		 analyses[:activity_gif]...)
end

# * Export
export analyse_WilsonCowan73_solution,

    find_local_maxima,

    plot_solution_surface
