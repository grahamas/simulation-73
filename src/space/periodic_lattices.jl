
abstract type AbstractPeriodicLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
@with_kw struct PeriodicLattice{T,N_ARR} <: AbstractLattice{T,N_ARR,N_ARR}
    extent::NTuple{N_ARR,T}
    n_points::NTuple{N_ARR,Int}
end
difference(p_lattice::PeriodicLattice, edge) = abs_difference_periodic(edge, p_lattice.extent)
Base.size(p_lattice::PeriodicLattice) = p_lattice.n_points

const Circle{T} = PeriodicLattice{T,1}
const Torus{T} = PeriodicLattice{T,2}
