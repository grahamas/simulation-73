module Meshes

using Markdown
using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update!

abstract type Space{T} <: Parameter{T} end

# * Segment
@with_kw struct Segment{DistT} <: Space{DistT}
    extent::DistT
    n_points::Int
end

function calculate(segment::Segment{T}) where T
    return LinRange{T}(-(segment.extent/2), (segment.extent/2), segment.n_points)
end

struct CalculatedSegment{DistT<:Number} <: CalculatedParam{Segment{DistT}}
    segment::Segment{DistT}
    value::LinRange{DistT}
    CalculatedSegment{DistT}(segment::Segment{DistT}) where DistT = new(segment, calculate(segment))
end

function Calculated(segment::Segment{T}) where T
    CalculatedSegment{T}(segment)
end

function update!(cs::CalculatedSegment, segment::Segment)
    if cs.segment == segment
        return false
    else
        cs.segment = segment
        cs.value = calculate(segment)
        return true
    end
end


import Base: step, zeros, length, size
step(cs::CalculatedSegment) = step(cs.value)
length(cs::CalculatedSegment) = length(cs.value)

size(s::Segment) = (s.n_points,)
size(cs::CalculatedSegment) = size(cs.value)
zeros(seg::Segment{T}) where T = zeros(T,seg.n_points)

export step, zeros, length

# * Distance matrix

@with_kw struct DistanceMatrix{T} <: Parameter{T}
    calculated_segment::CalculatedSegment{T}
end

function DistanceMatrix(seg::Segment{T}) where T
    DistanceMatrix{T}(Calculated(seg))
end

@doc doc"""
This matrix contains values such that the $j^{th}$ column of the $i^{th}$ row
contains the distance between locations $i$ and $j$ in the 1D space dimension provided.
"""
function distance_matrix(dm::DistanceMatrix{T}) where {T <: Real}
    # aka Hankel, but that method isn't working in SpecialMatrices
    xs = dm.calculated_segment.value
    distance_mx = zeros(T, length(xs), length(xs))
    for i in range(1, length=length(xs))
        distance_mx[:, i] = abs.(xs .- xs[i])
    end
    return distance_mx'
end

struct CalculatedDistanceMatrix{T} <: CalculatedParam{DistanceMatrix{T}}
    distance_matrix::DistanceMatrix{T}
    value::Matrix{T}
    CalculatedDistanceMatrix{T}(dm::DistanceMatrix{T}) where T = new(dm, distance_matrix(dm))
end
step(cdm::CalculatedDistanceMatrix) = step(cdm.distance_matrix.calculated_segment)

function CalculatedDistanceMatrix(segment::S) where {T,S <: Segment{T}}
    CalculatedDistanceMatrix{T}(DistanceMatrix(segment))
end

function CalculatedDistanceMatrix(calculated_segment::S) where {T,S <: CalculatedSegment{T}}
    CalculatedDistanceMatrix{T}(calculated_segment)
end

function Calculated(dm::DistanceMatrix{T}) where T
    CalculatedDistanceMatrix{T}(dm)
end


function update!(cdm::CalculatedDistanceMatrix, segment::Segment)
    if update!(cdm.calc_segment, segment)
        cdm.value = distance_matrix(cdm.calc_segment)
        return true
    else
        return false
    end
end

export DistanceMatrix, CalculatedDistanceMatrix

export Space, Segment, Calculated, update!

end
