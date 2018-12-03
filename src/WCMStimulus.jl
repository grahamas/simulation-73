module WCMStimulus

using Parameters

using Meshes
using CalculatedParameters
import CalculatedParameters: Calculated, update!

abstract type Stimulus{T} <: Parameter{T} end

#---------- CompoundStimulus ------------#

@with_kw struct AddedStimuli{T} <: Stimulus{T}
    stimuli::Array{Stimulus{T}} # Sadly can't specify all must be subtype of Stimulus
end

stim_param(stim::Stimulus{T}) where T = T

function AddedStimuli(stims...)
    T = stim_param(stims[1])
    stims = Union{map(typeof,stims)...}[stims...]
    AddedStimuli{T}(stims)
end

function add(arr::Array{Stimulus{T}}) where T
    AddedStimuli{T}(arr)
end

function add(arr::Array)
    pop_stims = zip(arr...)
    T = stim_param(collect(pop_stims)[1][1])
    return AddedStimuli{T}[AddedStimuli(stims...) for stims in pop_stims]
end

function calculate(as::AddedStimuli, space::Segment)
    calculated_stimuli = map((stim) -> Calculated(stim, space), as.stimuli)
    (args) -> sum(map((calculated_stimulus) -> get_value(calculated_stimulus)(args...), calculated_stimuli))
end

mutable struct CalculatedAddedStimuli{T} <: CalculatedParam{AddedStimuli{T}}
    stimulus::AddedStimuli{T}
    space::Segment{T}
    value::Function
    CalculatedAddedStimuli{TT}(s::AddedStimuli{TT}, space::Segment{TT}) where TT = new(s, space, calculate(s,space))
end

function Calculated(as::AddedStimuli{T}, space::Segment{T}) where T
    CalculatedAddedStimuli{T}(as, space)
end

export AddedStimuli, CalculatedAddedStimuli, add

# ------------- GaussianNoiseStimulus ----------- #

function gaussian_noise(space, mean, sd) # assumes signal power is 0db
    return mean .+ sd .* randn(size(space))
end

struct GaussianNoiseStimulus{T} <: Stimulus{T}
    mean::T
    sd::T
end

function GaussianNoiseStimulus{T}(; SNR::T=0.0, mean::T=0.0) where T
    sd = sqrt(1/10 ^ (SNR / 10))
    GaussianNoiseStimulus{T}(mean, sd)
end

function calculate(wns::GaussianNoiseStimulus, space::Segment)
    (t) -> gaussian_noise(space, wns.mean, wns.sd) # Not actually time dependent
end

mutable struct CalculatedGaussianNoiseStimulus{T} <: CalculatedParam{GaussianNoiseStimulus{T}}
    stimulus::GaussianNoiseStimulus{T}
    space::Segment{T}
    value::Function
    CalculatedGaussianNoiseStimulus{TT}(s::GaussianNoiseStimulus{TT}, space::Segment{TT}) where TT = new(s, space, calculate(s,space))
end

function Calculated(wns::GaussianNoiseStimulus{T}, space::Segment{T}) where T
    CalculatedGaussianNoiseStimulus{T}(wns, space)
end

export GaussianNoiseStimulus, CalculatedGaussianNoiseStimulus

# ----------- SharpBumpStimulus ------------ #

struct SharpBumpStimulus{T} <: Stimulus{T}
    width::T
    strength::T
    window::Tuple{T,T}
end

function SharpBumpStimulus{T}(; strength=nothing, width=nothing,
        duration=nothing, window=nothing) where T
    if window == nothing
        return SharpBumpStimulus{T}(width, strength, (zero(T), duration))
    else
        @assert duration == nothing
        return SharpBumpStimulus{T}(width, strength, window)
    end
end

function SharpBumpStimulus(p)
    T = typeof(p[:(Stimulus.width)])
    SharpBumpStimulus{T}(p[:(Stimulus.width)], p[:(Stimulus.strength)], p[:(Stimulus.window)])
end

function calculate(sbs::SharpBumpStimulus, space::Segment)
    sharp_bump_factory(Calculated(space), sbs.width, sbs.strength, sbs.window)
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
        csbs.value = calculate(sbs, csbs.space)
        return true
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
    frame = zero(mesh_coords)
    half_width = width / 2      # using truncated division
    xs = mesh_coords[:,1]   # Assumes all pops have same mesh_coords
    start_dx = findfirst(xs .>= mid_point - half_width)
    stop_dx = findlast(xs .<= mid_point + half_width)
    frame[start_dx:stop_dx,:] .= strength
    return frame
end
"""
The "sharp bump" is the usual theoretical impulse: Binary in both time and
space. On, then off.
"""
function sharp_bump_factory(segment::CalculatedParam{Segment{DistT}}, width, strength, (onset, offset)) where {DistT <: Real}
        # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_sharp_bump_frame(segment, width, strength)
    off_frame = zero(on_frame)
    return (t) -> (onset <= t < offset) ? on_frame : off_frame
end

end