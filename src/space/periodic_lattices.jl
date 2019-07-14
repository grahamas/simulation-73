
abstract type AbstractPeriodicLattice{T,D} <: AbstractLattice{T,D} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
@with_kw struct PeriodicLattice{T,N} <: AbstractLattice{T,N}
    extent::NTuple{N,T}
    n_points::NTuple{N,Int}
end
distance_metric(p_lattice::PeriodicLattice, edge) = euclidean_metric_periodic(edge, p_lattice.extent)
Base.size(p_lattice::PeriodicLattice) = p_lattice.n_points

const Circle{T} = PeriodicLattice{T,1}
const Torus{T} = PeriodicLattice{T,2}
