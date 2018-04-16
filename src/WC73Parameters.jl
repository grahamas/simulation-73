module WC73Parameters

using Parameters
import JSON

push!(LOAD_PATH, ".")
using SimulationTypes

# * Parameter object
# ** Definition
@with_kw struct WilsonCowan73Params{InteractionType, ParamType}
  # Explict fields in parameter file
  # May also be given as LaTeX command (e.g. alpha for α)
    α::ParamType     # Weight on homeostatic term
    β::ParamType     # Weight on nonlinear term
    τ::ParamType     # Time constant
    a::ParamType     # Sigmoid steepness
    θ::ParamType     # Sigmoid translation
    r::ParamType     # Refractory period multiplier
  # Other fields in parameter file include
  # :time => {[:N], :extent}
  # :space => {:N, :extent}
  # :stimulus => {:weight, :duration, :strength}
  # :connectivity => {:amplitudes, :spreads}
  # Constructed fields
    W::InteractionType    # Tensor interaction multiplier
    stimulus_fn::Function
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

    stimulus_params = expand_params(mesh, pop!(p, :stimulus))
    connectivity_params = expand_params(mesh, pop!(p, :connectivity))
    p = expand_params(mesh, p)

    p[:mesh] = mesh
    p[:stimulus_fn] = make_stimulus_fn(mesh; stimulus_params...)
    p[:W] = sholl_connectivity(mesh, connectivity_params[:amplitudes],
                               connectivity_params[:spreads])

    return WilsonCowan73Params(; p...)
end

function expand_params(mesh::AbstractMesh, dct::T) where T <: Dict
    for (k,v) in dct
        if v isa PopulationParam
            dct[k] = expand_param(mesh, v)
        end
    end
    return dct
end

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
    if a[1] isa Array && a[1][1] isa Number # eltype gives Any, for some reason
        return InteractionParam(vcat([convert_py(arr) for arr in a]...))
    elseif a[1] isa Dict
        return convert_py.(a)
    elseif a[1] isa Number
        return PopulationParam(convert_py.(vcat(a...))) # Python arrays are rows...
    else
        error("Unsupported parse input array of eltype $(typeof(a[1]))")
    end
end

convert_py(val::String) = val

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
    convert_pykey(k::String) = (convert_pykey ∘ Symbol)(k)

    return Dict(convert_pykey(k) => convert_py(v) for (k,v) in d)
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

# ** Loading function

function load_WilsonCowan73_parameters(json_filename::String, modifications=nothing)
    # Parse JSON with keys as symbols.
    param_dct = (convert_py ∘ JSON.parsefile)(json_filename)
    return deep_merge(param_dct, modifications)
end
# ** Export and END
export load_WilsonCowan73_parameters, WilsonCowan73Params
end
