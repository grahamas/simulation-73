abstract type AbstractSpace{T,D} <: AbstractParameter{T} end

discrete_segment(extent::T, n_points::Int) where {T <: Number} = LinRange{T}(-(extent/2), (extent/2), n_points)
@calculated_type(struct Circle{T} <: AbstractSpace{T,1}
    extent::T
    n_points::Int
end, function calculate()
    discrete_segment(extent, n_points)
end
)
distance_metric(circ::Circle, a, b) = min(mod(a-b, circ.extent), mod(b-a, circ.extent))


@calculated_type(struct Segment{T} <: AbstractSpace{T,1}
    extent::T
    n_points::Int
end, function calculate()
    discrete_segment(extent, n_points)
end
)
distance_metric(::Segment, a, b) = abs(a - b)

@calculated_type(struct Pops{P,T,D,S} <: AbstractSpace{T,D}
    space::S
end, function calculate()
    calculate(space)#repeat(calculate(space), outer=(ones(Int,D)...,P))
end
)

@generated function one_pop(calc_pops::CalculatedType{<:Pops{P,T,D,S}}) where {P,T,D,S}
    colons = [:(:) for i in 1:D]
    quote
        calc_pops.value[$(colons...),1]
    end
end

Pops{n_pops}(space::S) where {T,D,n_pops,S <: AbstractSpace{T,D}} = Pops{n_pops,T,D,S}(space)
distance_metric(ps::Pops, a, b) = distance_metric(ps.space, a, b)

@calculated_type(struct Grid{T} <: AbstractSpace{T,2}
    extent::Tuple{T,T}
    n_points::Tuple{Int,Int}
end, function calculate()
    zip(discrete_segment.(extent, n_points))
end
)

function get_edges(locations::CalculatedType{<:AbstractSpace{T}}) where {T}
    ((locations.value[i], locations.value[j]) for (i,j) in Iterators.product(CartesianIndices(locations.value),CartesianIndices(locations.value)))
end

#region Segment

get_space_origin_idx(line::AbstractSpace{T,1}) where T = CartesianIndex(round(Int, line.n_points / 2), 1)
get_space_origin_idx(pops::Pops) = get_space_origin_idx(pops.space)

import Base: step, zeros, length, size, ndims
step(line::AbstractSpace{T,1}) where T = line.extent / (line.n_points - 1)
length(line::AbstractSpace{T,1}) where T = line.n_points

step(ps::Pops) = step(ps.space)
length(ps::Pops) = length(ps.space)

step(cs::CalculatedType{<:AbstractSpace}) = step(cs.source)
length(cs::CalculatedType{<:AbstractSpace}) = length(cs.source)

size(cs::CalculatedType{<:AbstractSpace}) = size(cs.value)

size(line::AbstractSpace{T,1}) where T = (length(line),)
size(ps::Pops{P}) where P = (size(ps.space)...,P)

zeros(space::AbstractSpace{T}) where {T} = zeros(T,size(space)...)

ndims(space::AbstractSpace) = length(size(space))
#endregion

#
# Base.zero(tup::Type{T}) where {T <: Tuple} = T(map(zero, t.parameters))
# @calculated_type(struct DistanceMatrix{COORD_T,CS} <: AbstractParameter{COORD_T}
#     calculated_space::CS
# end,
# begin
#     xs = calculated_space.value[:,1]
#     distance_mx = zeros(COORD_T, length(xs), length(xs))
#     for i in range(1, length=length(xs))
#         distance_mx[:, i] .= distance_metric.(Ref(calculated_space.source), xs, xs[i])
#     end
#     return distance_mx'
# end,
# Array{T,2}
# )
#
# step(cdm::CalculatedDistanceMatrix) = step(cdm.source.calculated_space)
#
# function DistanceMatrix{COORD_T}(cs::CS) where {COORD_T,CS}
#     DistanceMatrix{COORD_T,CS}(cs)
# end
# function DistanceMatrix(space::S) where {T, S <: AbstractSpace{T}}
#     DistanceMatrix{T}(Calculated(space))
# end
# function CalculatedDistanceMatrix(space::S) where {T,S <: AbstractSpace{T}}
#     CalculatedDistanceMatrix{T}(DistanceMatrix(space))
# end
