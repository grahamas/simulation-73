module Meshes

using Markdown
using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update!

abstract type Space{T,N} <: Parameter{T} end
abstract type PopSpace{T,N,P} <: Space{T,N} end

# * Segment
@with_kw struct PopSegment{DistT,P} <: PopSpace{DistT,1,P}
    extent::DistT
    n_points::Int
end

function calculate(segment::PopSegment{T,P}) where {T,P}
    one_pop = LinRange{T}(-(segment.extent/2), (segment.extent/2), segment.n_points)
    return repeat(one_pop, outer=(1,P))
end

struct CalculatedPopSegment{DistT<:Number,P} <: CalculatedParam{PopSegment{DistT,P}}
    segment::PopSegment{DistT}
    value::Array{DistT,2}
    CalculatedPopSegment{DistT,P}(segment::PopSegment{DistT,P}) where {DistT,P} = new(segment, calculate(segment))
end

function Calculated(segment::PopSegment{T,P}) where {T,P}
    CalculatedPopSegment{T,P}(segment)
end

function update!(cs::CalculatedPopSegment, segment::PopSegment)
    if cs.segment == segment
        return false
    else
        cs.segment = segment
        cs.value = calculate(segment)
        return true
    end
end

get_origin(seg::PopSegment) = CartesianIndex(round(Int, seg.n_points / 2), 1)

import Base: step, zeros, length, size, ndims
step(seg::PopSegment) = seg.extent / (seg.n_points - 1)
step(cs::CalculatedPopSegment) = step(cs.segment)
length(cs::CalculatedPopSegment) = length(cs.value)

size(s::PopSegment{T,P}) where {T,P} = (s.n_points,P)
size(cs::CalculatedPopSegment) = size(cs.value)
zeros(seg::PopSegment{T,P}) where {T,P} = zeros(T,seg.n_points,P)

ndims(space::PopSpace{T,N}) where {T,N} = N + 1

export step, zeros, length, get_origin

# * Distance matrix

@with_kw struct DistanceMatrix{T} <: Parameter{T}
    calculated_segment::CalculatedPopSegment{T}
end

function DistanceMatrix(seg::PopSegment{T}) where T
    DistanceMatrix{T}(Calculated(seg))
end

@doc doc"""
This matrix contains values such that the $j^{th}$ column of the $i^{th}$ row
contains the distance between locations $i$ and $j$ in the 1D space dimension provided.
"""
function distance_matrix(dm::DistanceMatrix{T}) where {T <: Real}
    # aka Hankel, but that method isn't working in SpecialMatrices
    xs = dm.calculated_segment.value[:,1]
    distance_mx = zeros(T, length(xs), length(xs))
    for i in range(1, length=length(xs))
        distance_mx[:, i] .= abs.(xs .- xs[i])
    end
    return distance_mx'
end

struct CalculatedDistanceMatrix{T} <: CalculatedParam{DistanceMatrix{T}}
    distance_matrix::DistanceMatrix{T}
    value::Matrix{T}
    CalculatedDistanceMatrix{T}(dm::DistanceMatrix{T}) where T = new(dm, distance_matrix(dm))
end
step(cdm::CalculatedDistanceMatrix) = step(cdm.distance_matrix.calculated_segment)

function CalculatedDistanceMatrix(segment::S) where {T,S <: PopSegment{T}}
    CalculatedDistanceMatrix{T}(DistanceMatrix(segment))
end

function CalculatedDistanceMatrix(calculated_segment::S) where {T,S <: CalculatedPopSegment{T}}
    CalculatedDistanceMatrix{T}(calculated_segment)
end

function Calculated(dm::DistanceMatrix{T}) where T
    CalculatedDistanceMatrix{T}(dm)
end


function update!(cdm::CalculatedDistanceMatrix, segment::PopSegment)
    if update!(cdm.calc_segment, segment)
        cdm.value = distance_matrix(cdm.calc_segment)
        return true
    else
        return false
    end
end

export DistanceMatrix, CalculatedDistanceMatrix

export Space, PopSegment, PopSpace, PopSegment, Calculated, update!

end
