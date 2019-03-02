
# Subsampling space can definitely be necessary
# Subsampling time may not be, maybe rather use DESolution methods.

# Supposing adaptive time, what do I want to do?
#	- Restrict analysis resolution
#	- Regularize analysis interval (!!!)

abstract type AbstractSubsampler{DT,W} end

@with_kw struct Subsampler{DT,W} <: AbstractSubsampler{DT,W}
	Δ::DT = nothing
	window::W = nothing
end

@with_kw struct IndexInfo{D}
	Δ::D
	origin_idx::Int # TODO: Make CartesianIndex
end

struct StrideToEnd
	stride::Int
	start::Int
	StrideToEnd(stride, start=1) = new(stride, start)
end

import Base: to_indices, _maybetail, @_inline_meta, tail, getindex
to_indices(A, inds, I::Tuple{StrideToEnd, Vararg{Any}})=
	(@_inline_meta; (I[1].start:I[1].stride:inds[1][end], to_indices(A, _maybetail(inds), tail(I))...))
import Base: getindex
getindex(A, S::StrideToEnd) = getindex(A, to_indices(A, (S,))...)
#getindex(A, ind::StrideToEnd) = A[ind.start:ind.stride:end]

function scalar_to_idx_window(scalar_window::Tuple{T,T}, arr::AbstractArray{T,1}) where T
	(findfirst((scalar_window[1] .< arr) .| (scalar_window[1] .≈ arr)), findlast((scalar_window[2] .> arr) .| (scalar_window[2] .≈ arr)))
end

function subsampling_Δidx(Δsubsampled::T, Δsource::T) where T
	Δidx = max(1, round(Int, Δsubsampled / Δsource))
	return Δidx
end

function subsampling_idxs(Δsource::T, origin_idx::Int, Δsubsampled::T, scalar_window::Tuple{T,T}) where T
	Δidx = subsampling_Δidx(Δsubsampled, Δsource)
	if scalar_window[1] == -Inf
		lower_idx = 1
	else
		lower_idx = round(Int, origin_idx + scalar_window[1] / Δsource)
	end
	if scalar_window[2] == Inf
		return StrideToEnd(Δidx,lower_idx)
	else
		upper_idx = round(Int, origin_idx + scalar_window[2] / Δsource)
		return lower_idx:Δidx:upper_idx
	end
end

function subsampling_idxs(info, subsampler::Subsampler{Nothing, Nothing})
	Colon()
end

function subsampling_idxs(info::IndexInfo{T}, subsampler::Subsampler{T,Tuple{T,T}}) where T
	subsampling_idxs(info.Δ, info.origin_idx, subsampler.Δ, subsampler.window)
end

function subsampling_idxs(info::IndexInfo{T}, subsampler::Subsampler{T,Nothing}) where {T <: Number}#Δsubsampled::T, scalar_window::Nothing) where T
	Δidx = subsampling_Δidx(subsampler.Δ, info.Δ)
	return StrideToEnd(Δidx)
end

function subsampling_idxs(info::IndexInfo{T}, subsampler::Subsampler{Nothing,Tuple{T,T}}) where T
	Δsubsampled = info.Δ
	subsampling_idxs(info.Δ, info.origin_idx, Δsubsampled, subsampler.window)
end

function subsampling_idxs(target::AbstractArray{T,1}, source::AbstractArray{T,1}) where T
	Δsource = source[2] - source[1]

	Δsubsampled = target[2] - target[1]
	scalar_window = (target[1], target[end])
	idx_window = scalar_to_idx_window(scalar_window, source)

	Δidx = subsampling_Δidx(Δsubsampled, Δsource)

	idx_window[1]:Δidx:idx_window[2]
end

function subsampling_idxs(target::AbstractArray{T,1}, source_info::IndexInfo) where T
	Δ = target[1] - target[2]
	target_subsampler = Subsampler(; window = (target[1], target[end]), Δ = target[2]-target[1])
	@assert all((target[2:end] - target[1:end-1]) .≈ target_subsampler.Δ)
	subsampling_idxs(source_info, target_subsampler)
end
