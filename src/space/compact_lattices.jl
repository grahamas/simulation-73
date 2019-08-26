
abstract type AbstractCompactLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A lattice of points, each an `N_ARR` coordinate, evenly distributed.

Can be constructed using default AbstractLattice constructor
"""
struct CompactLattice{T,N_ARR} <: AbstractCompactLattice{T,N_ARR,N_ARR}
    arr::Array{NTuple{N_ARR,T},N_ARR}
end
difference(lattice::AbstractCompactLattice, edge) = abs_difference(edge)


# Use default lattice plotting
