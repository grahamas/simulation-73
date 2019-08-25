
abstract type AbstractCompactLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
@with_kw struct CompactLattice{T,N_ARR} <: AbstractCompactLattice{T,N_ARR,N_ARR}
    extent::NTuple{N_ARR,T}
    n_points::NTuple{N_ARR,Int}
end
difference(lattice::AbstractCompactLattice, edge) = abs_difference(edge)


# Use default lattice plotting
