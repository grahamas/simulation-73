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
