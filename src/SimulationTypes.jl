using Parameters

# * Aliases
const PopulationParam{T} = RowVector{T, Array{T,1}} where T <: Real
const InteractionParam{T} = Array{T, 2} where T <: Real
const ExpandedParam{T} = Array{T, 2} where T <: Real
const ExpandedParamFlat{T} = Array{T, 1} where T <: Real
const InteractionTensor{T} = Array{T, 4} where T <: Real
const Interaction1DFlat{T} = Array{T, 2} where T <: Real
const SpaceState1D{T} = Array{T, 2} where T <: Real
const SpaceState1DFlat{T} = Array{T, 1} where T <: Real
const SpaceDim{DistT} = StepRangeLen{DistT} where T <: Real

export PopulationParam, InteractionParam, ExpandedParam
export ExpandedParamFlat, InteractionTensor, Interaction1DFlat
export SpaceState1D, SpaceState1DFlat, SpaceDim

# * Structs
# The structure used to input parameters.
# This maybe should be in a different module...?

@with_kw struct SpaceDimParameters
    extent::Int
    N::Int
end

@with_kw struct SpaceParameters
    dimensions::Array{SpaceDimParameters}
end

@with_kw struct ConnectivityParameters{T}
    amplitudes::InteractionParam{T}
    spreads::InteractionParam{T}
end

@with_kw struct StimulusParameters{ValueT<:Real, TimeT<:Real, DistT<:Real}
    name::String
    strength::ValueT
    duration::TimeT
    width::DistT
end

@with_kw struct WC73ModelParameters{ValueT<:Real, TimeT<:Real, DistT<:Real}
    space::SpaceParameters
    connectivity::ConnectivityParameters{ValueT}
    stimulus::StimulusParameters{ValueT, TimeT, DistT}
    β::PopulationParam{ValueT}
    τ::PopulationParam{ValueT}
    α::PopulationParam{ValueT}
    r::PopulationParam{ValueT}
    a::PopulationParam{ValueT}
    θ::PopulationParam{ValueT}
end

@with_kw struct SolverParameters{TimeT<:Real}
    T::TimeT
    dt::TimeT
end

@with_kw struct AnalysesParameters
    simulation_name::String
    root::String
    activity_gif::Dict
end

@with_kw struct WC73InputParameters{ValueT<:Real, TimeT<:Real, DistT<:Real}
    model::WC73ModelParameters{ValueT, TimeT, DistT}
    solver::SolverParameters{TimeT}
    analyses::AnalysesParameters
end
export SpaceDimParameters, SpaceParameters, ConnectivityParameters, StimulusParameters
export WC73ModelParameters, SolverParameters, AnalysesParameters, WC73InputParameters

# * Mesh
# Define a mesh type that standardizes interaction with the discretization of
# space (and populations, though those are inherently discrete, as we currently
# conceptualize them).

# ** Type Definition and Constructors
# All meshes are subtyped from AbstractMesh. SpaceMesh contains only discretized
# spatial dimensions. PopMesh contains a SpaceMesh, but also an integer indicating
# the number of colocalized populations (i.e. each spatial point contains members
# of each population). FlatMesh is merely a flattened representation of a PopMesh
# containing only one spatial dimension. Rather than concatenating populations
# along a "population dimension," the populations are concatenated along the
# single spatial dimension. This is useful so that the convolution can be
# implemented as a matrix multiplication, however I don't see how to extend
# it. I would not have implemented it, except that's how the preceding Python
# implementation worked, and I needed to have a direct comparison in order to
# debug.

abstract type AbstractMesh end
struct SpaceMesh{DistT} <: AbstractMesh
    dims::Array{SpaceDim{DistT}}
end
struct PopMesh{DistT} <: AbstractMesh
    space::SpaceMesh{DistT}
    n_pops::Integer
end
struct FlatMesh{DistT} <: AbstractMesh
    pop_mesh::PopMesh{DistT}
    FlatMesh{DistT}(mesh) where {DistT <: Real} = ndims(mesh) != 2 ? error("cannot flatten >1D mesh.") : new(mesh)
end

# Flatten and unflatten take a PopMesh to a FlatMesh and vice versa (only if the
# PopMesh has only 1D space).
flatten(mesh::PopMesh) = FlatMesh(mesh)
unflatten(mesh::FlatMesh) = mesh.pop_mesh

# FlatMesh has no outer constructor, as it uses the more descriptive "flatten."
function SpaceMesh(DistT::DataType, dim_dcts::Array{D}) where {D <: Dict}
    dims = Array{StepRangeLen}(length(dim_dcts))
    for (i, dim) in enumerate(dim_dcts)
        extent::DistT = dim[:extent]
        N::Integer = dim[:N]
        dims[i] = linspace(-(extent/2), (extent/2), N)
    end
    SpaceMesh{DistT}(dims)
end
function PopMesh{DistT}(dim_dcts::Array{<:Dict}, n_pops::Integer) where {DistT <: Real}
    PopMesh(SpaceMesh(dim_dcts),n_pops)
end

export SpaceMesh

# ** Methods
# Numerous functions operating on meshes, including size, ndims, true_ndims,
# coords, zeros, and expand_param.

import Base: size, ndims, zeros
function size(mesh::SpaceMesh)
    return length.(mesh.dims)
end
function size(mesh::PopMesh)
    return (size(mesh.space)..., mesh.n_pops)
end
function size(mesh::FlatMesh)
    return size(mesh.pop_mesh)[1] * mesh.pop_mesh.n_pops
end
function ndims(mesh::AbstractMesh)
    return length(size(mesh))
end
function true_ndims(mesh::AbstractMesh)
    return ndims(mesh)
end
# true_ndims returns the "real" structure of the mesh, i.e. unflattened.
function true_ndims(mesh::FlatMesh)
    return ndims(mesh.pop_mesh)
end
function coords(mesh::SpaceMesh)
    @assert ndims(mesh) == 1
    return mesh.dims[1]
end
function coords(mesh::PopMesh)
    @assert ndims(mesh) == 2
    return repeat(coords(mesh.space), outer=(1, mesh.n_pops))
end
function coords(mesh::FlatMesh)
    return repeat(coords(mesh.pop_mesh.space), outer=mesh.pop_mesh.n_pops)
end
function zeros(mesh::AbstractMesh)
    zeros(coords(mesh))
end
function expand_param(mesh::PopMesh, param::RowVector)::ExpandedParam
    space_dims = size(mesh)[1:end-1]
    return repeat(param, inner=(space_dims..., 1))
end
function expand_param(mesh::FlatMesh, param::RowVector)::ExpandedParamFlat
    return expand_param(mesh.pop_mesh, param)[:]
end

# ** Interface for applying functions
function apply(fn, mesh::SpaceMesh)
    return hcat([fn.(dim) for dim in mesh.dims]...)
end

function apply_with_time(fn, mesh::SpaceMesh, time)
    return hcat([fn.(dim, time) for dim in mesh.dims]...)
end

function apply_through_time(fn, mesh::SpaceMesh, time_len, dt)
    time_range = 0:dt:time_len
    output = Array{Float64,2}(size(mesh)..., length(time_range))
    for (i_time, time) in enumerate(time_range)
        output[:, i_time] = apply_with_time(fn, mesh, time)
    end
    return output
end

export apply, apply_through_time
