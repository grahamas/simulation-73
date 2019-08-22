
abstract type AbstractLattice{T,N_ARR,N_CONN} <: AbstractSpace{T,N_ARR,N_CONN} end
Base.step(space::AbstractLattice{T}) where T = space.extent ./ (space.n_points .- 1)

"""
    discrete_segment(extent, n_points)

Return an object containing `n_points` equidistant coordinates of a segment of length `extent` centered at 0. If you want 0 to be an element of the segment, make sure `n_points` is odd.

# Example
```jldoctest
julia> seg = discrete_segment(5.0, 7);

julia> length(seg) == 5
true

julia> seg[end] - seg[1] ≈ 5.0
true
```
"""
function discrete_segment(start::T, stop::T, n_points::Int) where {T <: Number}
    LinRange{T}(start, stop, n_points)
end
function discrete_segment(extent::T, n_points::Int) where {T <: Number}
    discrete_segment(-extent/2,extent/2,n_points)
end
"""
    coordinates(extent, n_points)

Return an object containing `n_points` equidistant coordinates along each dimension of a grid of length `extent` along each dimension, centered at (0,0,...,0).
"""
function discrete_lattice(extent::NTuple{N,T}, n_points::NTuple{N,Int}) where {N,T}
    Iterators.product(discrete_segment.(extent, n_points)...) |> collect
end
coordinates(lattice::AbstractLattice) = discrete_lattice(lattice.extent, lattice.n_points)
coordinate_axes(lattice::AbstractLattice) = (discrete_segment.(lattice.extent, lattice.n_points)...,)
"""
    origin_idx(lattice)

Return the coordinate of the zero point in `lattice`.

# Example
```jldoctest
julia> segment = CompactLattice(10.0, 11)
CompactLattice{Float64,1}(10.0, 11)

julia> seg_origin_idx = origin_idx(segment)
CartesianIndex(6,)

julia> collect(coordinates(segment))[seg_origin_idx] == 0.0
true

julia> grid = CompactLattice((10.0,50.0), (11, 13))
PeriodicLattice{Float64,2}((10.0, 50.0), (11, 13))

julia> grid_origin_idx = origin_idx(CompactLattice((10.0,50.0), (11, 13)))
CartesianIndex(6, 7)

julia> collect(coordinates(grid))[grid_origin_idx] == (0.0, 0.0)
true
```
"""
origin_idx(lattice::AbstractLattice) = CartesianIndex(round.(Int, size(lattice) ./ 2, RoundNearestTiesUp))

using Statistics
@recipe function f(lattice::AbstractLattice, values::AbstractArray{<:AbstractArray,1})
    standard_deviations = std.(values)
    means = mean.(values)
    up_stds = means .+ standard_deviations
    down_stds = means .- standard_deviations
    seriestype := :line
    @series begin
        linestyle := [:dot :dot :solid]
        (lattice, hcat(up_stds, down_stds, means))
    end
end
