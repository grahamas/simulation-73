
"""
    `AbstractSpace{T,N_ARR,N_CDT}`
        with distance-type `T`,
        storable in array of dimension `N_ARR`,
        and with point-coordinates of dimension `N_CDT`.
    """
abstract type AbstractSpace{T,N_ARR,N_CDT} <: AbstractParameter{T} end

"""
    coordinates(space::AbstractSpace)

Return an object in the shape of the space where each element is the coordinate of that element.
"""
coordinates(space::AbstractSpace) = error("undefined.")

"""
    distances(space)

Return the distances between every pair of points in `space`
"""
function differences(space::AbstractSpace{T}) where T
    edges = Iterators.product(coordinates(space), coordinates(space))
    return difference.(Ref(space), edges)
end
function differences(space::AbstractSpace{T,N_ARR,N_CDT},
                     reference_location::NTuple{N_CDT,T}) where {T,N_ARR,N_CDT}
    edges = ((coord, reference_location) for coord in coordinates(space))
    return difference.(Ref(space), edges)
end

# Extend Base methods to AbstractSpace types
import Base: step, zero, length, size, ndims

zero(::Type{NTuple{N,T}}) where {N,T} = NTuple{N,T}(zero(T) for _ in 1:N)
zero(space::AbstractSpace{T}) where {T} = zeros(T,size(space)...)

ndims(space::AbstractSpace) = length(size(space))

# Include file definitions of various spaces
include("space/metrics.jl")
include("space/abstract_lattices.jl")
include("space/compact_lattices.jl")
include("space/periodic_lattices.jl")
include("space/lattice_augmentations.jl")
