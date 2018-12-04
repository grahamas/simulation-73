module WCMConnectivity

using Meshes
using CalculatedParameters
import CalculatedParameters: Calculated, update!
using Parameters

# * Types

abstract type Connectivity{T} <: Parameter{T} end

@with_kw struct ShollConnectivity{T} <: Connectivity{T}
    amplitude::T
    spread::T
end

function ShollConnectivity(p)
    ShollConnectivity(p[:("Connectivity.amplitude")], p[:("Connectivity.spread")])
end

mutable struct CalculatedShollConnectivity{T} <: CalculatedParam{ShollConnectivity{T}}
    connectivity::ShollConnectivity{T}
    calc_dist_mx::CalculatedDistanceMatrix{T}
    value::Matrix{T}
    function CalculatedShollConnectivity{T}(c::ShollConnectivity{T},d::CalculatedDistanceMatrix{T}) where T
        new(c, d, sholl_matrix(c, d))
    end
end

# Uhhhh why is this function not just inside Calculated???
# ShollConnectivity is specified by the argument; no need for name
function CalculatedShollConnectivity(connectivity::ShollConnectivity{T}, segment::Segment{T}) where T
    calc_dist_mx = Calculated(DistanceMatrix(segment))
    return CalculatedShollConnectivity{T}(connectivity, calc_dist_mx)
end

function Calculated(connectivity::ShollConnectivity, segment::Segment)
    CalculatedShollConnectivity(connectivity, segment)
end

function update!(csc::CalculatedShollConnectivity, connectivity::ShollConnectivity)
    if csc.connectivity == connectivity
        return false
    else
        csc.connectivity = connectivity
        csc.value = sholl_matrix(connectivity, csc.calc_dist_mx)
        return true
    end
end

function update!(csc::CalculatedShollConnectivity, connectivity::ShollConnectivity, space::Segment)
    if update!(csc.calc_dist_mx, space)
        csc.connectivity = connectivity
        csc.value = sholl_matrix(csc.connectivity, csc.calc_dist_mx)
        return true
    else
        return update!(csc, connectivity)
    end
end


# * Sholl connectivity

function sholl_matrix(connectivity::ShollConnectivity, calc_dist_mx::CalculatedDistanceMatrix)
    A = connectivity.amplitude
    σ = connectivity.spread
    dist_mx = calc_dist_mx.value
    step_size = step(calc_dist_mx)
    return sholl_matrix(A, σ, dist_mx, step_size)
    # In comparing to Neuman, there is no division by 2 on the edges
    # but the edges are very small, so there's no difference
end

"""
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

export Connectivity, ShollConnectivity, CalculatedShollConnectivity, Calculated, update!

end
