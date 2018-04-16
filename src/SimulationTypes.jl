using Parameters

# * SimulationTypes

const NumType = Float64
const DistType = NumType
const TimeType = NumType
const PopulationParam = RowVector{NumType, Array{NumType,1}}
const InteractionParam = Array{NumType, 2}
const ExpandedParam = Array{NumType, 2}
const ExpandedParamFlat = Array{NumType,1}
const InteractionTensor = Array{NumType,4}
const Interaction1DFlat = Array{NumType,2}
const SpaceState1D = Array{NumType, 2}
const SpaceState1DFlat = Array{NumType,1}
const SpaceDim = StepRangeLen{NumType}

export NumType, DistType, TimeType, PopulationParam, InteractionParam, ExpandedParam
export ExpandedParamFlat, InteractionTensor, Interaction1DFlat
export SpaceState1D, SpaceState1DFlat, SpaceDim

# The structure used to input parameters.
# This maybe should be in a different module...?

@with_kw struct SpaceDimParameters
    extent::Int
    N::Int
end

@with_kw struct SpaceParameters
    dimensions::Array{SpaceDimParameters}
end

@with_kw struct ConnectivityParameters
    amplitudes::InteractionParam
    spreads::InteractionParam
end

@with_kw struct StimulusParameters
    name::String
    duration::TimeType
    strength::NumType
    width::DistType
end

@with_kw struct WC73ModelParameters
    space::SpaceParameters
    connectivity::ConnectivityParameters
    stimulus::StimulusParameters
    β::PopulationParam
    τ::PopulationParam
    α::PopulationParam
    r::PopulationParam
    a::PopulationParam
    θ::PopulationParam
end

@with_kw struct SolverParameters
    T::TimeType
    dt::TimeType
end

@with_kw struct AnalysesParameters
    simulation_name::String
    root::String
    activity_gif::Dict
end

@with_kw struct WC73InputParameters
    model::WC73ModelParameters
    solver::SolverParameters
    analyses::AnalysesParameters
end
export SpaceDimParameters, SpaceParameters, ConnectivityParameters, StimulusParameters
export WC73ModelParameters, SolverParameters, AnalysesParameters, WC73InputParameters


# Define a mesh type that standardizes interaction with the discretization of
# space (and populations, though those are inherently discrete, as we currently
# conceptualize them).

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
struct SpaceMesh <: AbstractMesh
    dims::Array{SpaceDim}
end
struct PopMesh <: AbstractMesh
    space::SpaceMesh
    n_pops::Integer
end
struct FlatMesh <: AbstractMesh
    pop_mesh::PopMesh
    FlatMesh(mesh) = ndims(mesh) != 2 ? error("cannot flatten >1D mesh.") : new(mesh)
end

# Flatten and unflatten take a PopMesh to a FlatMesh and vice versa (only if the
# PopMesh has only 1D space).
flatten(mesh::PopMesh) = FlatMesh(mesh)
unflatten(mesh::FlatMesh) = mesh.pop_mesh

# FlatMesh has no outer constructor, as it uses the more descriptive "flatten."
function SpaceMesh(dim_dcts::Array{T}) where T <: Dict
    dims = Array{StepRangeLen}(length(dim_dcts))
    for (i, dim) in enumerate(dim_dcts)
        extent::NumType = dim[:extent]
        N::Integer = dim[:N]
        dims[i] = linspace(-(extent/2), (extent/2), N)
    end
    SpaceMesh(dims)
end
function PopMesh(dim_dcts::Array{T}, n_pops::Integer) where T <: Dict
    PopMesh(SpaceMesh(dim_dcts),n_pops)
end

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
