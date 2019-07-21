
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
"""
    coordinates(extent, n_points)

Return an object containing `n_points` equidistant coordinates along each dimension of a grid of length `extent` along each dimension, centered at (0,0,...,0).
"""
@memoize function coordinates(lattice::AbstractLattice)
     Iterators.product(coordinate_axes(lattice)...)
 end
function coordinate_axes(lattice::AbstractLattice)
    all(lattice.n_points .% 2 .== 1) || @warn "n_points = $n_points is not odd, so the segment will not have the origin."
    (discrete_segment.((-).(lattice.extent ./ 2), lattice.extent ./ 2, lattice.n_points)...,)
end
"""
    origin_idx(lattice)

Return the coordinate of the zero point in `lattice`.

# Example
```jldoctest
julia> segment = Segment(10.0, 11)
Segment{Float64}(10.0, 11)

julia> seg_origin_idx = origin_idx(segment)
CartesianIndex(6,)

julia> collect(coordinates(segment))[seg_origin_idx] == 0.0
true

julia> grid = Grid((10.0,50.0), (11, 13))
Grid{Float64}((10.0, 50.0), (11, 13))

julia> grid_origin_idx = origin_idx(Grid((10.0,50.0), (11, 13)))
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
