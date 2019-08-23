
abstract type AbstractCompactLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
@with_kw struct CompactLattice{T,N_ARR} <: AbstractCompactLattice{T,N_ARR,N_ARR}
    extent::NTuple{N_ARR,T}
    n_points::NTuple{N_ARR,Int}
end
difference(lattice::AbstractCompactLattice, edge) = abs_difference(edge)

@recipe function f(lattice::CompactLattice, values)
    x := coordinate_axes(lattice)[1] |> collect
    y := values
    seriestype := :line
    if val_lim != nothing
        ylim := val_lim
    end
    ()
end

@recipe function f(lattice::CompactLattice{T,2}, values::Array{T,2}; val_lim=nothing) where T
    (x, y) = coordinate_axes(lattice) .|> collect
    seriestype := :heatmap
    if val_lim != nothing
        clim := val_lim
        zlim := val_lim
    end
    aspect_ratio := :equal
    (x,y,values)
end
