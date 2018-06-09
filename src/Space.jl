module Space

# * Distance matrix
doc"""
This matrix contains values such that the $j^{th}$ column of the $i^{th}$ row
contains the distance between locations $i$ and $j$ in the 1D space dimension provided.
"""
function distance_matrix(xs::CalculatedSegment{T}) where {T <: Real}
    # aka Hankel, but that method isn't working in SpecialMatrices
    distance_mx = zeros(T, length(xs), length(xs))
    for i in range(1, length(xs))
        distance_mx[:, i] = abs.(xs - xs[i])
    end
    return distance_mx'
end

mutable struct CalculatedDistanceMatrix{T<:Number} <: Calculated{Matrix{T}}
    calc_segment::CalculatedSegment
    value::Matrix{T}
    CalculatedDistanceMatrix{T}(calc_seg::CalculatedSegment{T}) = new(calc_seg, distance_matrix(calc_seg))
end
step(cdm::CalculatedDistanceMatrix) = step(cdm.calc_segment)

function CalculatedDistanceMatrix(segment::Segment)
    CalculatedDistanceMatrix(CalculatedSegment(segment))
end
step(cs::CalculatedSegment) = step(cs.mesh)

function update!(cdm::CalculatedDistanceMatrix, segment::Segment)
    if update!(cdm.calc_segment, segment)
        cdm.value = distance_matrix(cdm.calc_segment)
        return true
    else
        return false
    end
end

# * Segment
struct Segment{DistT<:Number}
    extent::DistT
    n_points::Int
end

function calculate_segment(segment::Segment)
    return linspace(-(segment.extent/2), (segment.extent/2), segment.n_points)
end

mutable struct CalculatedSegment{DistT<:Number} <: Calculated{Segment{DistT}}
    segment::Segment{DistT}
    value::StepRangeLen{DistT}
    CalculatedSegment{DistT}(segment) = new(segment, calculate_segment(segment))
end

function update!(cs::CalculatedSegment, segment::Segment)
    if cs.segment == segment
        return false
    else
        cs.segment = segment
        cs.value = calculate_segment(segment)
        return true
    end
end



end
