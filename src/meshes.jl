"AbstractSpace{T,D} with distance type T and dimension D"
abstract type AbstractSpace{T,D} <: AbstractParameter{T} end

"""
    coordinates(space::AbstractSpace)
    coordinates(calc_space::CalculatedType{<:AbstractSpace})

Return an object in the shape of the space where each element is the coordinate of that element.
"""
coordinates(space::AbstractSpace) = coordinates(calculate(space))
coordinates(calc_space::CalculatedType{<:AbstractSpace}) = calc_space.value

"""
    euclidean_metric(edge)

Return the distance between two points in euclidean space, given an edge between those points.
"""
euclidean_metric(edge::Tuple{T,T}) where T<:Number = abs(edge[1] - edge[2])
euclidean_metric(edge::Tuple{Tup,Tup}) where {T,N,Tup<:NTuple{N,T}} = abs.(edge[1] .- edge[2])
"""
    periodic_euclidean_metric(edge, period)

Return the distance between two points in euclidean space as in euclidean_metric, but let the space wrap with period.
"""
function periodic_euclidean_metric(edge::Tuple{T,T}, period::T) where T<:Number
    diff = euclidean_metric(edge)
    if diff > period / 2
        return period - diff
    else
        return diff
    end
end
function periodic_euclidean_metric(edge::Tuple{Tup,Tup}, periods::Tup) where {N,T,Tup<:NTuple{N,T}}
    diffs = euclidean_metric(edge)
    diffs = map(zip(diffs, periods)) do (diff, period)
        if diff > period / 2
            return period - diff
        else
            return diff
        end
    end
    return Tup(diffs)
end

"""
    discrete_segment(extent, n_points)

Return an object containing `n_points` equidistant coordinates of a segment of length `extent` centered at 0.
"""
discrete_segment(extent::T, n_points::Int) where {T <: Number} = LinRange{T}(
    -(extent/2),(extent/2), n_points)
"""
    discrete_grid(extent, n_points)

Return an object containing `n_points` equidistant coordinates along each dimension of a grid of length `extent` along each dimension, centered at (0,0,...,0).
"""
discrete_grid(extent::NTuple{N,T}, n_points::NTuple{N,Int}) where {N,T} = Iterators.product(
    discrete_segment.(extent, n_points)...)

"A simple Segment of length `extent` and `n_points`-many points"
@calculated_type(struct Segment{T} <: AbstractSpace{T,1}
    extent::T
    n_points::Int
end, function calculate()
    discrete_segment(extent, n_points)
end
)
distance_metric(segment::Segment, edge) = euclidean_metric(edge)

"A Circle with circumference `extent` represented by `n_points`-many points"
@calculated_type(struct Circle{T} <: AbstractSpace{T,1}
    extent::T
    n_points::Int
end, function calculate()
    discrete_segment(extent, n_points)
end
)
distance_metric(circle::Circle, edge) = periodic_euclidean_metric(edge, circle.extent)

"""
    Pops{P}(space)

Wrap a generic AbstractSpace so that it repeats P-many times
"""
@calculated_type(struct Pops{P,T,D,S} <: AbstractSpace{T,D}
    space::S
end, function calculate()
    calculate(space)#repeat(calculate(space), outer=(ones(Int,D)...,P))
end
)
Pops{n_pops}(space::S) where {T,D,n_pops,S <: AbstractSpace{T,D}} = Pops{n_pops,T,D,S}(space)
distance_metric(pops::Pops, edge) = distance_metric(pops.space, edge)

"""
    one_pop(calc_pops)

Return the coordinates for a single population of a multi-Pops space.
"""
function one_pop(calc_pops::CalculatedType{<:Pops{P,T,D,S}}) where {P,T,D,S}
    coordinates(calc_pops)
end

"""A square Grid of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension."""
@calculated_type(struct Grid{T} <: AbstractSpace{T,2}
    extent::Tuple{T,T}
    n_points::Tuple{Int,Int}
end, function calculate()
    discrete_grid(extent, n_points)
end
)
distance_metric(grid::Grid, edge) = euclidean_metric(edge)

"""A Torus shaped grid of points."""
@calculated_type(struct Torus{T} <: AbstractSpace{T,2}
    extent::Tuple{T,T}
    n_points::Tuple{Int,Int}
end, function calculate()
    discrete_grid(extent, n_points)
end
)
distance_metric(torus::Torus, edge) = periodic_euclidean_metric(edge, torus.extent)

"""
    get_distances(calc_space)

Return the distances between every pair of points in `calc_space`
"""
function get_distances(calc_space::CalculatedType{<:AbstractSpace{T}}) where T
    edges = Iterators.product(calc_space.value, calc_space.value)
    distances = distance_metric.(Ref(calc_space.source), edges)
end

"""
    get_space_origin_idx(space)

Return the coordinate of the zero point in `space`. Without respect to populations.

# Example
```jldoctest
julia> segment = Segment(10.0, 11)
Segment{Float64}(10.0, 11)

julia> origin_idx = Simulation73.get_space_origin_idx(segment)
CartesianIndex(6,)

julia> collect(Simulation73.calculate(segment))[origin_idx] == 0.0
true

julia> Simulation73.get_space_origin_idx(Pops{2}(segment))
CartesianIndex(6,)

julia> grid = Grid((10.0,50.0), (11, 13))
Grid{Float64}((10.0, 50.0), (11, 13))

julia> origin_idx = Simulation73.get_space_origin_idx(Grid((10.0,50.0), (11, 13)))
CartesianIndex(6, 7)

julia> collect(Simulation73.calculate(grid))[origin_idx] == (0.0, 0.0)
true
```
"""
get_space_origin_idx(space::AbstractSpace{T}) where T = CartesianIndex(round.(Int, space.n_points ./ 2, RoundNearestTiesUp))
get_space_origin_idx(pops::Pops) = get_space_origin_idx(pops.space)

# Extend Base methods to AbstractSpace types
import Base: step, zero, length, size, ndims
step(line::AbstractSpace{T,1}) where T = line.extent / (line.n_points - 1)
step(space::AbstractSpace{T}) where T = space.extent ./ (space.n_points .- 1)
length(line::AbstractSpace{T,1}) where T = line.n_points

step(ps::Pops) = step(ps.space)
length(ps::Pops) = length(ps.space)



size(space::AbstractSpace{T,1}) where T = (space.n_points,)
size(space::AbstractSpace{T}) where T = space.n_points
size(ps::Pops{P}) where P = (size(ps.space)...,P)

zero(::Type{NTuple{N,T}}) where {N,T} = NTuple{N,T}(zero(T) for i in 1:N)
zero(space::AbstractSpace{T}) where {T} = zeros(T,size(space)...)
one_pop_zero(pops::Pops) = zero(pops.space)
one_pop_zero(cp::CalculatedType{<:Pops}) = one_pop_zeros(cp.source)

step(cs::CalculatedType{<:AbstractSpace}) = step(cs.source)
length(cs::CalculatedType{<:AbstractSpace}) = length(cs.source)
size(cs::CalculatedType{<:AbstractSpace}) = size(cs.value)

ndims(space::AbstractSpace) = length(size(space))
