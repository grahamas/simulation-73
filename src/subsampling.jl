
# Subsampling space can definitely be necessary
# Subsampling time may not be, maybe rather use DESolution methods.

# Supposing adaptive time, specify "saveat" -- simple.
# Optional, allow specifying "tstops" which forces non-interpolation.


import Base: to_indices, _maybetail, @_inline_meta, tail, getindex

abstract type AbstractSubsampler{D} end
function subsample(obj, sub_arr::AbstractArray)
    for sub in sub_arr
        obj = subsample(obj, sub)
    end
    return obj
end
function subsample(lattice::LATTICE, sub::AbstractSubsampler{D}) where {T,D,LATTICE<:AbstractLattice{T,D}}
    idxs = coordinate_indices(lattice, sub)
    LATTICE(lattice.arr[idxs])
end
subsample(space::AbstractSpace, ::Nothing) = space
function getindex(lattice::LATTICE, sub::AbstractSubsampler{D}) where {T,D,LATTICE<:AbstractLattice{T,D}}
    idxs = coordinate_indices(lattice, sub)
    LATTICE(lattice.arr[idxs])
end

@with_kw struct IndexSubsampler{D} <: AbstractSubsampler{D}
	strides::NTuple{D,Int}
end
@with_kw struct RightCutFromValue{D,T} <: AbstractSubsampler{D}
	cut::NTuple{D,T}
end
export RightCutFromValue

struct RadialSlice <: AbstractSubsampler{2} end
"Subsample a space down to a new space"
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
# struct StrideToEnd{D}
#     stride::CartesianIndex{D}
# 	start::CartesianIndex{D}
# 	StrideToEnd(stride::CartesianIndex{D}, start::CartesianIndex{D}=CartesianIndex{D}(1)) where D = new{D}(stride, start)
# end
#end
to_indices(A, inds, I::Tuple{StrideToEnd, Vararg{Any}})=
	(@_inline_meta; (I[1].start:I[1].stride:inds[1][end], to_indices(A, _maybetail(inds), tail(I))...))
to_indices(A, inds, I::Tuple{NTuple{N,StrideToEnd}, Vararg{Any}}) where N = to_indices(A, inds, (I[1]..., _maybetail(I)...))
getindex(A::AbstractArray, S::StrideToEnd) = getindex(A, to_indices(A, (S,))...)

# function getindex(A::AbstractArray{T,N}, VW::ValueWindower{N,T}) where {T,N}
# 	value_window = VW.window
# 	(findfirst((value_window[1] .< arr) .| (value_window[1] .≈ arr)), findlast((value_window[2] .> arr) .| (value_window[2] .≈ arr)))
# end

## RadialSlice
"Get the coordinates of the subsample space within the larger space."
function coordinate_indices(space::AbstractLattice{T,2}, subsampler::RadialSlice) where T
	origin = origin_idx(space)
	[CartesianIndex(origin[1], x) for x in origin[2]:size(space,2)]
end

function coordinate_indices(lattice::AbstractLattice{T}, subsampler::RightCutFromValue{D,T}) where {T,D}
    axis_dxs = map(zip(subsampler.cut, coordinate_axes(lattice))) do (value, axis)
        findfirst(axis .>= value)
    end
    stride = 1
    cuts = Tuple(axis_dxs)
    return StrideToEnd.(stride, cuts)
end

function coordinate_indices(::Any, subsampler::IndexSubsampler{D}) where D
    return StrideToEnd.(subsampler.strides, 1)
end