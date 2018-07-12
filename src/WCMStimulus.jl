module WCMStimulus

using Parameters

using Meshes
using CalculatedParameters
import CalculatedParameters: Calculated, update!

abstract type Stimulus{T} <: Parameter{T} end

@with_kw struct SharpBumpStimulus{T} <: Stimulus{T}
    width::T
    strength::T
    duration::T
end

function SharpBumpStimulus(p)
    T = typeof(p[:(Stimulus.width)])
    SharpBumpStimulus{T}(p[:(Stimulus.width)], p[:(Stimulus.strength)], p[:(Stimulus.duration)])
end

function calculate(sbs::SharpBumpStimulus, space::Segment)
    sharp_bump_factory(Calculated(space), sbs.width, sbs.strength, sbs.duration)
end

mutable struct CalculatedSharpBumpStimulus{T} <: CalculatedParam{SharpBumpStimulus{T}}
    stimulus::SharpBumpStimulus{T}
    space::Segment{T}
    value::Function
    CalculatedSharpBumpStimulus{S}(s::SharpBumpStimulus{S}, space::Segment{S}) where S = new(s, space, calculate(s, space))
end

function Calculated(sbs::SharpBumpStimulus{T}, space::Segment{T}) where T
    CalculatedSharpBumpStimulus{T}(sbs, space)
end

function update!(csbs::CalculatedSharpBumpStimulus, sbs::SharpBumpStimulus)
    if csbs.stimulus == sbs
        return false
    else
        csbs.stimulus = sbs
        csbs.value = calculate(sbs, space)
    end
end
export Stimulus, SharpBumpStimulus, Calculated, update!

# * Stimulus functions

# ** Sharp bump
# TODO Understand these functions again.....
"Implementation of sharp_bump_frame used in sharp_bump_factory"
function make_sharp_bump_frame(segment::CalculatedParam{Segment{DistT}}, width::DistT, strength) where {DistT <: Real}
    make_sharp_bump_frame(segment.value, width, strength)
end

function make_sharp_bump_frame(mesh_coords::AbstractArray{DistT}, width::DistT, strength) where {DistT <: Real}
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
"""
The "sharp bump" is the usual theoretical impulse: Binary in both time and
space. On, then off.
"""
function sharp_bump_factory(segment::CalculatedParam{Segment{DistT}}, width, strength, duration) where {DistT <: Real}
        # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_sharp_bump_frame(segment, width, strength)
    off_frame = zeros(on_frame)
    return (t) -> (t <= duration) ? on_frame : off_frame
end

end