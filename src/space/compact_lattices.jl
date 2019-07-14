
abstract type AbstractCompactLattice{T,D} <: AbstractLattice{T,D} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
@with_kw struct CompactLattice{T,N} <: AbstractCompactLattice{T,N}
    extent::NTuple{N,T}
    n_points::NTuple{N,Int}
end
distance_metric(lattice::AbstractCompactLattice, edge) = euclidean_metric(edge)
Base.size(lattice::CompactLattice) = lattice.n_points

const Segment{T} = CompactLattice{T,1}
const Grid{T} = CompactLattice{T,2}
