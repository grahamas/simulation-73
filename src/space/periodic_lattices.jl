
abstract type AbstractPeriodicLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
@with_kw struct PeriodicLattice{T,N_ARR} <: AbstractPeriodicLattice{T,N_ARR,N_ARR}
    extent::NTuple{N_ARR,T}
    n_points::NTuple{N_ARR,Int}
end
coordinate_axes(p_lattice::PeriodicLattice) = (discrete_segment.(0.0, p_lattice.extent .- step(p_lattice), p_lattice.n_points)...,)
difference(p_lattice::PeriodicLattice, edge) = abs_difference_periodic(edge, p_lattice.extent)

@recipe function f(lattice::PeriodicLattice{T,1}, values) where T
    θ = coordinate_axes(lattice)[1] .* (2π / lattice.extent[1])
    #y = values .* sin.(θ)
    #x = values .* cos.(θ)
    seriestype := :scatter
    projection := :polar

    aspect_ratio := :equal

    (θ |> collect, r_values)
end

@recipe function f(lattice::PeriodicLattice{T,2}, values::Array{T,2}) where T
    (x, y) = coordinate_axes(lattice) .|> collect
    seriestype := :heatmap
    (x,y,values)
end
