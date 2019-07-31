
abstract type AbstractPeriodicLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
struct PeriodicLattice{T,N_ARR} <: AbstractLattice{T,N_ARR,N_ARR}
    extent::NTuple{N_ARR,T}
    n_points::NTuple{N_ARR,Int}
end
PeriodicLattice(; extent=nothing, n_points=nothing) = PeriodicLattice(extent, n_points)
difference(p_lattice::PeriodicLattice, edge) = abs_difference_periodic(edge, p_lattice.extent)
Base.size(p_lattice::PeriodicLattice) = p_lattice.n_points

@recipe function f(lattice::PeriodicLattice{T,1}, values) where T
    θ = coordinate_axes(lattice)[1] .* (2π / lattice.extent[1])
    #y = values .* sin.(θ)
    #x = values .* cos.(θ)
    seriestype := :scatter
    projection := :polar
    (collect(θ), values .+ 1.0)
end
