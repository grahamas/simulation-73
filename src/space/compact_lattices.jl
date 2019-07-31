
abstract type AbstractCompactLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
@with_kw struct CompactLattice{T,N_ARR} <: AbstractCompactLattice{T,N_ARR,N_ARR}
    extent::NTuple{N_ARR,T}
    n_points::NTuple{N_ARR,Int}
end
difference(lattice::AbstractCompactLattice, edge) = abs_difference(edge)
Base.size(lattice::CompactLattice) = lattice.n_points

@recipe function f(lattice::CompactLattice, values)
    x := coordinate_axes(lattice)[1] |> collect
    y := values
    seriestype := :line
    ()
end

@recipe function f(lattice::CompactLattice{T,2}, values::Array{T,2}) where T
    (x, y) = coordinate_axes(lattice) .|> collect
    seriestype := :heatmap
    (x,y,values)
end
