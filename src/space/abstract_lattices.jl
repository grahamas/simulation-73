
abstract type AbstractLattice{T,N_ARR,N_CDT} <: AbstractSpace{T,N_ARR,N_CDT} end
(t::Type{<:AbstractLattice})(n_points::Tuple, extent::Tuple) = t(discrete_lattice(extent, n_points))
(t::Type{<:AbstractLattice{T,1}})(n_points::Number,extent::Number) where T = t((n_points,),(extent,))
(t::Type{<:AbstractLattice})(; n_points, extent) = t(n_points, extent)

Base.step(space::AbstractLattice) = extent(space) ./ (size(space) .- 1)
Base.size(lattice::AbstractLattice) = size(lattice.arr)
Base.size(lattice::AbstractLattice, d::Int) = size(lattice.arr, d)
Base.zeros(lattice::AbstractLattice{T}) where T = zeros(T,size(lattice)...)
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
extent(lattice::AbstractLattice) = lattice.arr[end] .- lattice.arr[1]
coordinates(lattice::AbstractLattice) = discrete_lattice(extent(lattice), size(lattice))
coordinate_axes(lattice::AbstractLattice) = (discrete_segment.(extent(lattice), size(lattice))...,)


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


@recipe function f(lattice::AbstractLattice{T,1}, values; val_lim=nothing) where T
    x := coordinate_axes(lattice)[1] |> collect
    y := values
    seriestype := :line
    if val_lim != nothing
        ylim := val_lim
    end
    ()
end

@recipe function f(lattice::AbstractLattice{T,2}, values::Array{T,2}; val_lim=nothing) where T
    (x, y) = coordinate_axes(lattice) .|> collect
    seriestype := :heatmap
    if val_lim != nothing
        clim := val_lim
        zlim := val_lim
    end
    aspect_ratio := :equal
    (x,y,values)
end
