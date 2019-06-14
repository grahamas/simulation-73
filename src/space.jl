
"`AbstractSpace{T,D}` with distance-type `T` and dimension `D`"
abstract type AbstractSpace{T,D} <: AbstractParameter{T} end
abstract type AbstractLattice{T,D} <: AbstractSpace{T,D} end

"""
    coordinates(space::AbstractSpace)

Return an object in the shape of the space where each element is the coordinate of that element.
"""
coordinates(space::AbstractSpace) = error("undefined.")
get_space_origin_idx(space::AbstractSpace{T}) where T = CartesianIndex(round.(Int, space.n_points ./ 2, RoundNearestTiesUp))
