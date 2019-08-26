
# Subsampling space can definitely be necessary
# Subsampling time may not be, maybe rather use DESolution methods.

# Supposing adaptive time, specify "saveat" -- simple.
# Optional, allow specifying "tstops" which forces non-interpolation.

abstract type AbstractSubsampler{D} end
subsample(space::AbstractSpace, ::Nothing) = space

@with_kw struct IndexSubsampler{D} <: AbstractSubsampler{D}
	strides::NTuple{D,Int}
end
@with_kw struct ValueSubsampler{D,T} <: AbstractSubsampler{D}
	Δ::NTuple{D,T}
end
@with_kw struct ValueWindower{D,T} <: AbstractSubsampler{D}
	window::NTuple{D,T}
end
struct RadialSlice <: AbstractSubsampler{2} end
function subsample(lattice::AbstractLattice{T,2}, rs::RadialSlice) where T
	dxs = coordinate_indices(lattice, rs)
	slice_coords = lattice.arr[dxs]
	CompactLattice{T,1}([(coord[2],) for coord in slice_coords])
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

## RadialSlice
function coordinate_indices(space::AbstractLattice{T,2}, subsampler::RadialSlice) where T
	origin = origin_idx(space)
	[CartesianIndex(origin[1], x) for x in origin[2]:size(space,2)]
end
