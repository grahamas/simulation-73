
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

const Segment{T} = CompactLattice{T,1}
const Grid{T} = CompactLattice{T,2}

@recipe function f(lattice::Segment, values; val_lim = nothing)
    x := coordinate_axes(lattice)[1] |> collect
    y := values
    seriestype := :line
    if val_lim != nothing
        ylim := val_lim
    end
    ()
end

@recipe function f(lattice::Grid{T}, values::Array{T,2}; val_lim=nothing) where T
    (x, y) = coordinate_axes(lattice) .|> collect
    seriestype := :heatmap
    if val_lim != nothing
        clim := val_lim
        zlim := val_lim
    end
    (x,y,values)
end
