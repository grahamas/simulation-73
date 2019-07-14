
abstract type AbstractAugmentedLattice{T,D,L} <: AbstractLattice{T,D} end

struct RandomlyEmbeddedLattice{T,D,L<:AbstractLattice{T,D},E<:AbstractSpace{T},SUM_N} <: AbstractAugmentedLattice{T,D,L}
    lattice::L
    embedded_lattice::E
    coordinates::Array{NTuple{SUM_N,T},D}
end
function RandomlyEmbeddedLattice(; lattice::L, embedded_lattice::E) where {T,D,L<:AbstractLattice{T,D},E<:AbstractSpace{T}}
    embedded_coordinates = embed_randomly(lattice, embedded_lattice)
    RandomlyEmbeddedLattice(lattice, embedded_lattice, embedded_coordinates)
end
function embed_randomly(lattice, embedded_lattice)
    [(lattice_coord..., sample(embedded_lattice)...) for lattice_coord in coordinates(lattice)]
end
function sample(lattice::AbstractLattice)
    (rand(length(lattice.extent)...) .* lattice.extent) .- (lattice.extent ./ 2)
end
