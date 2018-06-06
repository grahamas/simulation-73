using Parameters
using JSON

# * Parameter object
# ** Definition
@with_kw struct WilsonCowan73Params{InteractionType, ParamType}
  # Explict fields in parameter file
  # May also be given as LaTeX command (e.g. alpha for α)
    α::ParamType     # Weight on homeostatic term
    β::ParamType     # Weight on nonlinear term
    τ::ParamType     # Time constant
    r::ParamType     # Refractory period multiplier
  # Other fields in parameter file include
  # :time => {[:N], :extent}
  # :space => {:N, :extent}
  # :stimulus => {:weight, :duration, :strength}
  # :connectivity => {:amplitudes, :spreads}
  # Constructed fields
    W::InteractionType    # Tensor interaction multiplier
    stimulus_fn::Function
    nonlinearity_fn::Function
    mesh::AbstractMesh
end
# ** Constructor and helper
function WilsonCowan73Params(p)
    p = deepcopy(p) # to prevent mutation
    npops = length(p[:r])

    space_dims = pop!(p, :space)
    @assert length(space_dims) == 1      # Currently only supports 1D
    mesh = PopMesh(space_dims, npops)
    if ndims(mesh) == 2
        mesh = flatten(mesh)
    end
    @assert mesh isa FlatMesh

    #stimulus_params = expand_param(mesh, pop!(p, :stimulus))
    #connectivity_params = expand_param(mesh, pop!(p, :connectivity))
    stimulus_params = pop!(p, :stimulus)
    connectivity_params = pop!(p, :connectivity)
    nonlinearity_params = pop!(p, :nonlinearity) #DO NOT EXPAND_PARAM
    p = expand_param(mesh, p)

    p[:mesh] = mesh
    p[:nonlinearity_fn] = make_nonlinearity_fn(mesh; nonlinearity_params...)
    p[:stimulus_fn] = make_stimulus_fn(mesh; stimulus_params...)
    p[:W] = sholl_connectivity(mesh, connectivity_params[:amplitudes],
                               connectivity_params[:spreads])
    return WilsonCowan73Params(; p...)
end

function expand_param(mesh::AbstractMesh, dct::D) where D <: Dict
    new_dct = Dict{Symbol,Any}()
    for (k,v) in dct
        new_dct[k] = expand_param(mesh, v)
    end
    return new_dct
end

function expand_param(mesh::AbstractMesh, v::Number)
    return v
end

function expand_param(mesh::AbstractMesh, s::AbstractString)
    return s
end

# ** Export
export WilsonCowan73Params
# * Load Parameters
#=
Because I originally wrote this in Python, the parameter files are JSON. (In the
process of moving to fully Julia parameters).
=#
# ** Conversion from Python to Julia
function convert_py(val::Number)
    float(val)
end

function convert_py(a::T) where T <: Array
    if a[1] isa Array# && a[1][1] isa Number # eltype gives Any, for some reason
        return vcat([convert_py(arr) for arr in a]...)
    elseif a[1] isa Dict
        return convert_py.(a)
    elseif a[1] isa Number
        return RowVector(convert_py.(vcat(a...))) # Python arrays are rows...
    elseif a[1] isa AbstractString
        return convert_py.(a)
    else
        error("Unsupported parse input array of eltype $(typeof(a[1]))")
    end
end

convert_py(val::String) = val

function finish(val::Array{T,1}) where {T <: Number}
    return RowVector(val)
end

function finish(val)
    return val
end

function convert_pyval(val)
    return (finish ∘ convert_py)(val)
end

function convert_py(d::T) where T <: Dict
    # TODO: Find package that does this...
    unicode_dct = Dict(:alpha=>:α, :beta=>:β, :tau=>:τ, :theta=>:θ)
    function convert_pykey(k_sym::Symbol)
        if k_sym in keys(unicode_dct)
            return unicode_dct[k_sym]
        else
            return k_sym
        end
    end
    function convert_pykey(k::String)
        (convert_pykey ∘ Symbol)(k)
    end
    return Dict(convert_pykey(k) => convert_pyval(v) for (k,v) in d)
end

function make_sure_dict(arr::Array)
    Dict(arr)
end

function make_sure_dict(dct::Dict)
    dct
end

# ** Merge dictionaries
function deep_merge(dct1, dct2::D) where D <: Dict
    new_dct = deepcopy(dct1)
    for k in keys(dct2)
        if k in keys(dct1)
            new_dct[k] = deep_merge(dct1[k], dct2[k])
        else
            new_dct[k] = dct2[k]
        end
    end
    return new_dct
end
function deep_merge(el1, el2)
    return el2
end
function deep_merge(el1, void::Void)
    return el1
end

# ** Standardizations
function standardize_parameters!(param_dct::Dict{Symbol}, json_filename::String)
    standardizations = [
        function check_simulation_name!(dct)
            if length(dct[:output][:simulation_name]) == 0
                dct[:output][:simulation_name] = splitext(basename(json_filename))[1]
            end
        end
    ]
    standardizations .|> (fn) -> fn(param_dct)
end

# ** Loading function

function load_WilsonCowan73_parameters(json_filename::String, modifications::Dict=nothing)
    # Parse JSON with keys as symbols.
    param_dct::Dict = (convert_py ∘ make_sure_dict ∘ JSON.parsefile)(json_filename)
    param_dct = deep_merge(param_dct, modifications)
    standardize_parameters!(param_dct, json_filename)
    return param_dct
end

# ** Export
export load_WilsonCowan73_parameters
# * Stimulus functions
doc"""
    make_stimulus_fn(mesh; name, stimulus_args...)

A factory taking the domain (`mesh`) and `name` of a stimulus and returning the function
defined to be associated with that name mapped over the given domain.
"""
function make_stimulus_fn(mesh; name=nothing, stimulus_args...)
    stimulus_factories = Dict(
        "smooth_bump" => smooth_bump_factory,
        "sharp_bump" => sharp_bump_factory,
        "two_sharp_bumps" => two_sharp_bumps_factory
    )
    return stimulus_factories[name](mesh; stimulus_args...)
end

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
function sharp_bump_factory(mesh; width=nothing, strength=nothing, duration=nothing)
        # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_sharp_bump_frame(mesh, width, strength)
    off_frame = zeros(on_frame)
    return (t) -> (t <= duration) ? on_frame : off_frame
end
