
abstract type AbstractPeriodicLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
struct PeriodicLattice{T,N_ARR} <: AbstractLattice{T,N_ARR,N_ARR}
    extent::NTuple{N_ARR,T}
    n_points::NTuple{N_ARR,Int}
end
PeriodicLattice(; extent=nothing, n_points=nothing) = PeriodicLattice(extent, n_points)
coordinate_axes(p_lattice::PeriodicLattice) = (discrete_segment.(0.0, p_lattice.extent .- step(p_lattice), p_lattice.n_points)...,)
difference(p_lattice::PeriodicLattice, edge) = abs_difference_periodic(edge, p_lattice.extent)
Base.size(p_lattice::PeriodicLattice) = p_lattice.n_points
origin_idx(p_lattice::PeriodicLattice{T,N_ARR}) where{T,N_ARR} = CartesianIndex(zero(T) for _ in 1:N)

const Circle{T} = PeriodicLattice{T,1}
const Torus{T} = PeriodicLattice{T,2}

@recipe function f(lattice::Circle{T}, values::AbstractArray{T}; val_lim = extrema(values)) where T
    θ = coordinate_axes(lattice)[1] .* (2π / lattice.extent[1])

    r_min, r_max = abs.(val_lim)
    if r_min < r_max
        r_axis = max(2 * r_min, r_max)
    else
        r_axis = 2 * r_min
    end

    lim = r_axis + r_max * 1.2
    ylim := (-lim, lim)
    xlim := (-lim, lim)

    r_values = hcat(values, zeros(size(θ)...)) .+ r_axis

    legend := :none
    projection := :polar

    (θ |> collect, r_values)
end
