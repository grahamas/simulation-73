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

    stimulus_params = expand_param(mesh, pop!(p, :stimulus))
    #connectivity_params = expand_param(mesh, pop!(p, :connectivity))
    connectivity_params = pop!(p, :connectivity)
    nonlinearity_params = expand_param(mesh, pop!(p, :nonlinearity))
    p = expand_param(mesh, p)

    p[:mesh] = mesh
    p[:nonlinearity_fn] = make_nonlinearity_fn(; nonlinearity_params...)
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
    if a[1] isa Array && a[1][1] isa Number # eltype gives Any, for some reason
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

# ** Standardizations
function standardize_parameters!(param_dct::Dict{Symbol}, json_filename::String)
    standardizations = [
        function check_simulation_name!(dct)
            if length(dct[:analyses][:simulation_name]) == 0
                dct[:analyses][:simulation_name] = splitext(basename(json_filename))[1]
            end
        end
    ]
    standardizations .|> (fn) -> fn(param_dct)
end

# ** Loading function

function load_WilsonCowan73_parameters(json_filename::String, modifications=nothing)
    # Parse JSON with keys as symbols.
    param_dct = (convert_py ∘ JSON.parsefile)(json_filename)
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
        "sharp_bump" => sharp_bump_factory
    )
    return stimulus_factories[name](mesh; stimulus_args...)
end

# ** Smooth bump
"Implementation of smooth_bump_frame used in smooth_bump_factory."
function make_smooth_bump_frame(mesh_coords::Array{DistT}, width::DistT, strength::ValueT, steepness::ValueT) where {ValueT <: Real, DistT <: Real}
    sig_diff_fn = make_sigmoid_diff_fn(; a=steepness, θ=-width/2, width=width)
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
function make_sharp_bump_frame(mesh::PopMesh{ValueT}, width::DistT, strength::ValueT) where {ValueT <: Real, DistT <: Real}
    mesh_coords = coords(mesh)
    frame = zeros(mesh_coords)
    mid_point = 0     # half length, half width
    half_width = width / 2      # using truncated division
    xs = mesh_coords[:,1]   # Assumes all pops have same mesh_coords
    start_dx = find(xs .>= mid_point - half_width)[1]
    stop_dx = find(xs .<= mid_point + half_width)[end]
    frame[start_dx:stop_dx,:] = strength
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

# * Connectivity functions

doc"""
This matrix contains values such that the $j^{th}$ column of the $i^{th}$ row
contains the distance between locations $i$ and $j$ in the 1D space dimension provided.
"""
function distance_matrix(xs::SpaceDim{DistT}) where {DistT <: Real}
    # aka Hankel, but that method isn't working in SpecialMatrices
    distance_mx = zeros(DistT, length(xs), length(xs))
    for i in range(1, length(xs))
        distance_mx[:, i] = abs.(xs - xs[i])
    end
    return distance_mx'
end
doc"""
We use an exponential connectivity function, inspired both by Sholl's
experimental work, and by certain theoretical considerations.

The interaction between two populations is entirely characterized by this
function and its two parameters: the amplitude (weight) and the spread
(σ). The spatial step size is also a factor, but as a computational concern
rather than a fundamental one.
"""
function sholl_matrix(amplitude::ValueT, spread::ValueT,
                      dist_mx::Interaction1DFlat{ValueT}, step_size::ValueT) where {ValueT <: Real}
    conn_mx = @. amplitude * step_size * exp(
        -abs(dist_mx / spread)
    ) / (2 * spread)
    return conn_mx
end
doc"""
This calculates a matrix of Sholl's exponential decay for each pair of
populations, thus describing all pairwise interactions. The result is a tensor
describing the effect of the source population at one location on the target
population in another location (indexed: `[tgt_loc, tgt_pop, src_loc,
src_pop]`). This works for arbitrarily many populations (untested) but only for
1D space.
"""
function sholl_connectivity(mesh::PopMesh{DistT}, W::InteractionParam{ValueT},
			    Σ::InteractionParam{ValueT})::InteractionTensor{ValueT} where {ValueT <: Real, DistT <: Real}
    xs = mesh.space.dims[1]
    N_x = length(xs)
    N_pop = size(W)[1]
    conn_tn = zeros(N_x, N_pop, N_x, N_pop)
    for tgt_pop in range(1,N_pop)
	for src_pop in range(1,N_pop)
	    conn_tn[:, tgt_pop, :, src_pop] .= sholl_matrix(W[tgt_pop, src_pop],
			  Σ[tgt_pop, src_pop], distance_matrix(xs), step(xs))
	end
    end
    return conn_tn
end
doc"""
In the two population case, flattening the tensor and using matrix
multiplication is 3x faster. This function provides exactly that.
"""
function sholl_connectivity(mesh::FlatMesh{ValueT}, args...) where {ValueT <: Real}
    # Why didn't I provide an unflattened mesh in the first place?
    sholl_connectivity(unflatten(mesh), args...) |> flatten_sholl
end
function flatten_sholl(tensor)::Interaction1DFlat
    N_x, N_p = size(tensor)[1:2]
    @assert N_p < N_x
    @assert size(tensor) == (N_x, N_p, N_x, N_p)
    flat = zeros(eltype(tensor), N_x*N_p, N_x*N_p)
    for i in 1:N_p
        for j in 1:N_p
            flat[(1:N_x)+((i-1)*N_x), (1:N_x)+((j-1)*N_x)] = tensor[:,i,:,j]
        end
    end
    return flat
end
# * Nonlinearity functions

function make_nonlinearity_fn(; name=error("Missing arg"), args=error("Missing arg"))
    nonlinearity_factories = Dict(
        "sigmoid" => make_sigmoid_fn,
        "sech2" => make_sech2_fn,
        "sigmoid_diff" => make_sigmoid_diff_fn

    )
    nonlinearity_factory = nonlinearity_factories[name]
    return nonlinearity_factory(; args...)
end

function rectify(x)
    return max(0,x)
end

# ** Sech2 functions

function make_sech2_fn(; a=error("Missing arg"), θ=error("Missing arg"))
    return (x) -> max.(0,sech2_fn(x, a, θ))
end

function sech2_fn(x, a, θ)
    return @. 1 - tanh(a * (x - θ))^2
end

# ** Sigmoid functions
function make_sigmoid_fn(; a=error("Missing arg"), θ=error("Missing arg"))
    return (x) -> sigmoid_fn(x, a, θ)
end

doc"""
The sigmoid function is defined
```math
\begin{align}
\mathcal{S}(x) = \frac{1}{1 + \exp(-a(x - \theta))}
\end{align}
```
where $a$ describes the slope's steepness and $\theta$ describes translation of the slope's center away from zero.

This is "simple" because in practice we use the rectified sigmoid.
"""
function simple_sigmoid_fn(x, a, theta)
    return @. (1 / (1 + exp(-a * (x - theta))))
end

doc"""
A rectified version of `simple_sigmoid_fn`.

In practice, we use rectified sigmoid functions because firing rates cannot be negative.

TODO: Rename to rectified_sigmoid_fn.
"""
function sigmoid_fn(x, a, theta)
    return max.(0, simple_sigmoid_fn(x, a, theta) .- simple_sigmoid_fn(0, a, theta))
end

# ** Difference of sigmoids functions
function make_sigmoid_diff_fn(; a=nothing, θ=nothing, width=nothing)
    unscaled(y) = sigmoid_diff_fn(y, a, θ, width)  # Peak is not always 1
    range = (θ-(1.0 ./ a)):0.001:(θ+(1.0 ./ a)+width)
    println("range $range")
    maxes = maximum(unscaled.(range), 1)
    return (x) -> unscaled(x) ./ maxes
end

function sigmoid_diff_fn(input, a, θ, width)
    return max.(0,sigmoid_fn(input, a, θ) - sigmoid_fn(input, a, θ + width))
end
