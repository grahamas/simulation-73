module WCMTarget

using WC73
using Targets
using Parameters
import DifferentialEquations: DESolution

function penalize_increasing(timeseries)
	difference = timeseries[2:end] .- timeseries[1:end-1]
	increases = difference[difference .> 0]
	return sum(increases)
end

function penalize_decreasing(timeseries)
	difference = timeseries[2:end] .- timeseries[1:end-1]
	decreases = difference[difference .<= 0]
	return sum(abs.(decreases))
end

function long_halflife(peak_vals)
	beginning = peak_vals[1]
	posthalf = find(peak_vals .< beginning/2)
	if length(posthalf) < 1
		return Inf
	end
	return 1.0 / posthalf[1]
end

function must_travel(peak_vals, peak_dxs::Array{Int},
	small_number=0.01, small_dx=27)
	normal_peaks = peak_vals .> small_number
	normal_peak_dxs = peak_dxs[normal_peaks]
	if sum(normal_peak_dxs .> small_dx) < 1
		return Inf
	end
	return 0
end

function must_not_increase(timeseries)
	difference = timeseries[2:end] .- timeseries[1:end-1]
	increases = difference[(difference .> 0)]
	if sum(increases) > 0
		return Inf
	end
	return sum(increases)
end

function must_decrease(timeseries)
	difference = timeseries[2:end] .- timeseries[1:end-1]
	decreases = difference[difference .< 0]
	if sum(decreases) == 0
		return Inf
	end
	return sum(decreases)
end

function must_be_pulse(timeseries)
	if timeseries[1] ≈ 0
		first_nonzero = find(timeseries .> 0)[1]
		timeseries = timeseries[first_nonzero:end]
	end
	if (length(timeseries) == 0) || (sum(timeseries .≈ 0) == 0)
		return Inf
	end
	return 0
end

@with_kw struct DecayingTraveling{T} <: LossFunction{T}
	timepoints::AbstractArray{T}
	target_pop::Int
	space_start::T
end

@with_kw struct Traveling{T} <: LossFunction{T}
	timepoints::AbstractArray{T}
	target_pop::Int
	space_start::T
end

function (p::DecayingTraveling{T})(soln::DESolution, calc_space) where {T}
	time_vec = p.timepoints
	pop = p.target_pop
	timeseries = soln(time_vec)[calc_space .> p.space_start,pop,:]
	peak_vals = Array{T,1}(length(time_vec))
	peak_dxs = Array{Int,1}(length(time_vec))
	for time_dx in 1:length(time_vec)
		peak_vals[time_dx], peak_dxs[time_dx] = findmax(timeseries[:,time_dx])
	end

	return (penalize_decreasing(peak_dxs)
			+ must_travel(peak_vals, peak_dxs)
			+ must_not_increase(peak_vals)
			+ must_decrease(peak_vals)
			+ must_be_pulse(timeseries[1,:]))
end

function (p::Traveling{T})(soln::DESolution, calc_space) where {T}
	time_vec = p.timepoints
	pop = p.target_pop
	timeseries = soln(time_vec)[calc_space .> p.space_start,pop,:]
	peak_vals = Array{T,1}(length(time_vec))
	peak_dxs = Array{Int,1}(length(time_vec))
	for time_dx in 1:length(time_vec)
		peak_vals[time_dx], peak_dxs[time_dx] = findmax(timeseries[:,time_dx])
	end

	return (penalize_decreasing(peak_dxs)
			+ must_travel(peak_vals, peak_dxs)
			+ must_not_increase(peak_vals)
			+ must_be_pulse(timeseries[1,:]))
end

@with_kw struct DecayingWaveFactory{T} <: TargetFactory{T}
    timepoints::AbstractArray
    decay::T
    speed::T
    target_pop::Int
end

function (F::DecayingWaveFactory)(model::WCMSpatial1D)
    segment = Calculated(model.space).value
    wave_fn(t) = @. exp(F.decay * t) * sech(segment - F.speed * t)
    l_wave_frames = wave_fn.(F.timepoints)
    target_pop = cat(3, l_wave_frames...)
    return cat(2, target_pop, zeros(target_pop))
end

export DecayingWaveFactory, DecayingTraveling

end