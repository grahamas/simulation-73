module Stimulus

using Parameters

using ..Space
using ..CalculatedParameter

type Stimulus{T} <: Parameter{T} end

@with_kw struct SharpBumpStimulus{T} <: Stimulus{T}
    width::T
    strength::T
    duration::T
end

function SharpBumpStimulus(p)
    SharpBumpStimulus(p[:(Stimulus.width)], p[:(Stimulus.strength)], p[:(Stimulus.duration)])
end

function calculate(sbs::SharpBumpStimulus, space::Space)
    sharp_bump_factory(calculate(space), sbs.width, sbs.strength, sbs.duration)
end

function Calculated(sbs::SharpBumpStimulus, space::Space)
    CalculatedSharpBumpStimulus(sbs, space, calculate(sbs, space))
end

mutable struct CalculatedSharpBumpStimulus{T} <: Calculated{SharpBumpStimulus}
    stimulus::SharpBumpStimulus{T}
    space::Segment
    value::Function
    CalculatedSharpBumpStimulus{T}(s::SharpBumpStimulus{T}) = new(s, calculate(s))
end

function update!(csbs::CalculatedSharpBumpStimulus, sbs::SharpBumpStimulus)
    if csbs.stimulus == sbs
        return false
    else
        csbs.stimulus = sbs
        csbs.value = calculate(stimulus, space)
    end
end


# * Stimulus functions

# ** Collision

function two_sharp_bumps_factory(mesh; width=nothing, strength=nothing, duration=nothing)
    @assert all([width, strength, duration] .!= nothing)
    on_frame = make_two_sharp_bumps_frame(mesh, width, strength)
    off_frame = zeros(on_frame)
    return (t) -> (t <= duration) ? on_frame : off_frame
end

function make_two_sharp_bumps_frame(mesh::PopMesh, args...)
    make_two_sharp_bumps_frame(coords(mesh), args...)
end

function make_two_sharp_bumps_frame(mesh::FlatMesh, args...)
    structured_frame = make_two_sharp_bumps_frame(mesh.pop_mesh, args...)
    return structured_frame[:]
end

function make_two_sharp_bumps_frame(mesh_coords::Array{DistT}, width::DistT, strength::ValueT) where {ValueT <: Real, DistT <: Real}
    frame = zeros(mesh_coords)
    mid_dx = floor(Int, size(mesh_coords,1)/2)
    frame[1:mid_dx,:] = make_sharp_bump_frame(mesh_coords[1:mid_dx,:], width, strength)
    frame[mid_dx+1:end,:] = make_sharp_bump_frame(mesh_coords[mid_dx+1:end,:], width, strength)
    return frame
end

# ** Smooth bump
"Implementation of smooth_bump_frame used in smooth_bump_factory."
function make_smooth_bump_frame(mesh_coords::Array{DistT}, width::DistT, strength::ValueT, steepness::ValueT) where {ValueT <: Real, DistT <: Real}
    mid_dx = floor(Int, length(mesh_coords) / 2)
    mid_value = mesh_coords(mid_dx)
    sig_diff_fn = make_neg_domain_sigmoid_diff_fn(steepness, mid_value-width/2, width)
    normed_bump_frame = sig_diff_fn.(mesh_coords)
    return strength .* normed_bump_frame
    #@. strength * (simple_sigmoid_fn(mesh_coords, steepness, -width/2) - simple_sigmoid_fn(mesh_coords, steepness, width/2))
end

"""
The smooth bump is a smooth approximation of the sharp impulse defined
elsewhere. It is smooth in both time and space. It is constructed essentially
from three sigmoids: Two coplanar in space, and one orthogonal to those in
time. The two in space describe a bump: up one sigmoid, then down a negative
sigmoid. The one in time describes the decay of that bump.

This stimulus has the advantages of being 1) differentiable, and 2) more
realistic. The differentiabiilty may be useful for the automatic solvers that
Julia has, which can try to automatically differentiate the mutation function
in order to improve the solving.
"""
function smooth_bump_factory(mesh::AbstractMesh;
                             width=nothing, strength=nothing, duration=nothing,
                             steepness=nothing)
    # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_smooth_bump_frame(coords(mesh), width, strength, steepness)
    return (t) -> @. on_frame * (1 - simple_sigmoid_fn(t, steepness, duration))
end

# ** Sharp bump
# TODO Understand these functions again.....
"Implementation of sharp_bump_frame used in sharp_bump_factory"
function make_sharp_bump_frame(mesh::PopMesh{ValueT}, width::DistT, strength) where {ValueT <: Real, DistT <: Real}
    make_sharp_bump_frame(coords(mesh), width, strength)
end
function make_sharp_bump_frame(segment::Calculated{Segment{DistT}}, width::DistT, strength) where {DistT <: Number}
    make_sharp_bump_frame(segment.value, width, strength)
end

function make_sharp_bump_frame(mesh_coords::Array{DistT}, width::DistT, strength::Union{ValueT,PopulationParam{ValueT}}) where {ValueT <: Real, DistT <: Real}
    mid_dx = floor(Int, size(mesh_coords, 1) / 2)
    mid_point = mesh_coords[mid_dx,1]
    frame = zeros(mesh_coords)
    half_width = width / 2      # using truncated division
    xs = mesh_coords[:,1]   # Assumes all pops have same mesh_coords
    start_dx = find(xs .>= mid_point - half_width)[1]
    stop_dx = find(xs .<= mid_point + half_width)[end]
    frame[start_dx:stop_dx,:] .= strength
    return frame
end
function make_sharp_bump_frame(mesh::FlatMesh, args...)
    structured_frame = make_sharp_bump_frame(mesh.pop_mesh, args...)
    flat_frame = structured_frame[:] # Works because FlatMesh must have 1D PopMesh
    return flat_frame
end
"""
The "sharp bump" is the usual theoretical impulse: Binary in both time and
space. On, then off.
"""
function sharp_bump_factory(segment::Calculated{Segment}, width, strength, duration)
        # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_sharp_bump_frame(segment, width, strength)
    off_frame = zeros(on_frame)
    return (t) -> (t <= duration) ? on_frame : off_frame
end
