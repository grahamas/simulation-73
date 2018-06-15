module Connectivity

using ..Space
using ..CalculatedParameter

# * Types

abstract type Connectivity end

@with_kw struct ShollConnectivity{T<:Number} <: Connectivity
    amplitude::T
    spread::T
end

function ShollConnectivity(p)
    ShollConnectivity(p[:("Connectivity.amplitude")], p[:("Connectivity.spread")])
end

mutable struct CalculatedShollConnectivity{T<:Number} <: Calculated{ShollConnectivity}
    connectivity::ShollConnectivity{T}
    calc_dist_mx::CalculatedDistanceMatrix{T}
    value::Matrix{T}
    CalculatedShollConnectivity{T}(c::ShollConnectivity{T},d::CalculatedDistanceMatrix{T}) = new(c, d, make_sholl_mx(c, d))
end

function CalculatedShollConnectivity(connectivity::ShollConnectivity, segment::Segment)
    calc_dist_mx = CalculatedDistanceMatrix(segment)
    return CalculatedShollConnectivity(connectivity, calc_dist_mx)
end

function Calculated(connectivity::ShollConnectivity, segment::Segment)
    CalculatedShollConnectivity(connectivity, segment)
end

function update!(csc::CalculatedShollConnectivity, connectivity::ShollConnectivity)
    if csc.connectivity == connectivity
        return false
    else
        csc.connectivity = connectivity
        csc.value = make_sholl_mx(connectivity, csc.calc_dist_mx)
        return true
    end
end

function update!(csc::CalculatedShollConnectivity, connectivity::ShollConnectivity, space::Space)
    if update!(csc.calc_dist_mx, space)
        csc.connectivity = connectivity
        csc.value = make_sholl_mx(csc.connectivity, csc.calc_dist_mx)
        return true
    else
        return update!(csc, connectivity)
    end
end


# * Sholl connectivity

function sholl_matrix(connectivity::ShollConnectivity, calc_dist_mx::CalculatedDistMatrix)
    A = connectivity.amplitude
    σ = connectivity.spread
    dist_mx = calc_dist_mx.value
    step_size = calc_dist_mx.step
    return sholl_matrix(A, σ, dist_mx, step_size)
end

doc"""
We use an exponential connectivity function, inspired both by Sholl's
experimental work, and by certain theoretical considerations.

The interaction between two populations is characterized by this
function and its two parameters: the amplitude (weight) and the spread
(σ). The spatial step size is also a factor, but as a computational concern
rather than a fundamental one.
"""
function sholl_matrix(amplitude::ValueT, spread::ValueT,
                      dist_mx::Array{ValueT,2}, step_size::ValueT) where {ValueT <: Real}
    @. amplitude * step_size * exp(
        -abs(dist_mx / spread)
    ) / (2 * spread)
end

end
