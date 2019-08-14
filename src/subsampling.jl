
# Subsampling space can definitely be necessary
# Subsampling time may not be, maybe rather use DESolution methods.

# Supposing adaptive time, specify "saveat" -- simple.
# Optional, allow specifying "tstops" which forces non-interpolation.

abstract type AbstractSubsampler{D} end

@with_kw struct IndexSubsampler{D} <: AbstractSubsampler{D}
	strides::NTuple{D,Int}
end
@with_kw struct ValueSubsampler{D,T} <: AbstractSubsampler{D}
	Δ::NTuple{D,T}
end
@with_kw struct ValueWindower{D,T} <: AbstractSubsampler{D}
	window::NTuple{D,T}
end

"IndexInfo specifies the index of 0 and the stride (`Δ`) between indices."
@with_kw struct IndexInfo{D,N}
	Δ::NTuple{N,D}
	origin_idx::CartesianIndex{N}
end

"""
	StrideToEnd(stride, start=1)

StrideToEnd is a custom index type that acts like `start:stride:end`, to circumvent the fact that you can't put `start:stride:end` into a variable.
"""
struct StrideToEnd
    stride::Int
	start::Int
	StrideToEnd(stride::Int, start::Int=1) = new(stride, start)
end
import Base: to_indices, _maybetail, @_inline_meta, tail, getindex
to_indices(A, inds, I::Tuple{StrideToEnd, Vararg{Any}})=
	(@_inline_meta; (I[1].start:I[1].stride:inds[1][end], to_indices(A, _maybetail(inds), tail(I))...))
getindex(A::AbstractArray, S::StrideToEnd) = getindex(A, to_indices(A, (S,))...)

getindex(A::AbstractArray{T,N}, IS::IndexSubsampler{N}) where {T,N} = A[[StrideToEnd(stride) for stride in IS.strides]...]
function getindex(A::AbstractArray{T,N}, VW::ValueWindower{N,T}) where {T,N}
	value_window = VW.window
	(findfirst((value_window[1] .< arr) .| (value_window[1] .≈ arr)), findlast((value_window[2] .> arr) .| (value_window[2] .≈ arr)))
end
function getindex(A::AbstractArray{T,N}, VS::ValueSubsampler{N,T}) where {T,N}


# function subsampling_Δidx(Δsubsampled::T, Δsource::T) where T
# 	Δidx = max(1, round(Int, Δsubsampled / Δsource))
# 	return Δidx
# end
#
#
# function subsampling_idxs(Δsource::T, origin_idx::Int, Δsubsampled::T, scalar_window::Tuple{T,T}) where {N,T<:Number}
# 	Δidx = subsampling_Δidx(Δsubsampled, Δsource)
# 	lower_idx = max(1, round(Int, origin_idx + (scalar_window[1] / Δsource)))
# 	if scalar_window[2] == Inf
# 		return StrideToEnd(Δidx,lower_idx)
# 	else
# 		upper_idx = round(Int, origin_idx + scalar_window[2] / Δsource)
# 		return lower_idx:Δidx:upper_idx
# 	end
# end
#
# function subsampling_idxs(Δsource::T, origin_idx::CartesianIndex{N}, Δsubsampled::T, scalar_window::Tuple{T,T}) where {N,T<:NTuple{N}}
# 	[subsampling_idxs(args...) for args in zip(Δsource, Tuple(origin_idx), Δsubsampled, zip(scalar_window...))]
# end
#
# function subsampling_idxs(info, subsampler::Subsampler{Nothing, Nothing})
# 	Colon()
# end
#
# function subsampling_idxs(info::IndexInfo{T,N}, subsampler::Subsampler{T,Tuple{T,T}}) where {N, T<:NTuple{N}}
# 	subsampling_idxs(info.Δ, info.origin_idx, subsampler.Δ, subsampler.window)
# end
#
# function subsampling_idxs(info::IndexInfo{T}, subsampler::Subsampler{T,Nothing}) where {T <: Number}#Δsubsampled::T, scalar_window::Nothing) where T
# 	Δidx = subsampling_Δidx(subsampler.Δ, info.Δ)
# 	return StrideToEnd(Δidx)
# end
#
# function subsampling_idxs(info::IndexInfo{T}, subsampler::Subsampler{Nothing,Tuple{T,T}}) where T
# 	Δsubsampled = info.Δ
# 	subsampling_idxs(info.Δ, info.origin_idx, Δsubsampled, subsampler.window)
# end

# FIXME generalize to N-D
# function subsampling_idxs(target::AbstractArray{T,1}, source::AbstractArray{T,1}) where T
# 	Δsource = (source[2] - source[1],)
#
# 	Δsubsampled = (target[2] - target[1],)
# 	scalar_window = (target[1], target[end])
# 	idx_window = scalar_to_idx_window(scalar_window, source)
#
# 	Δidx = subsampling_Δidx(Δsubsampled, Δsource)
#
# 	idx_window[1]:Δidx:idx_window[2]
# end
#
# function subsampling_idxs(target::AbstractArray{T,1}, source_info::IndexInfo) where T
# 	Δ = target[1] - target[2]
# 	target_subsampler = Subsampler(; window = (target[1], target[end]), Δ = target[2]-target[1])
# 	@assert all((target[2:end] - target[1:end-1]) .≈ target_subsampler.Δ)
# 	subsampling_idxs(source_info, target_subsampler)
# end
