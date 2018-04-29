module WC73
using Parameters
import JSON

using TensorOperations

using DifferentialEquations

ENV["GKSwstype"] = "100" # For headless plotting (on server)
using Plots; gr()

import Base.Dates
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
function load_WilsonCowan73_parameters(json_filename::String, modifications::Union{Dict,Void})
    # Parse JSON with keys as symbols.
    param_dct = (convert_py ∘ JSON.parsefile)(json_filename)
    return deep_merge(param_dct, modifications)
end

const NumType = Float64
const DistType = NumType
const TimeType = NumType
const PopulationParam = RowVector{NumType}
const InteractionParam = Array{NumType, 2}
const ExpandedParam = Array{NumType, 2}
const ExpandedParamFlat = Array{NumType,1}
const InteractionTensor = Array{NumType,4}
const Interaction1DFlat = Array{NumType,2}
const SpaceState1D = Array{NumType, 2}
const SpaceState1DFlat = Array{NumType,1}
const SpaceDim = StepRangeLen{NumType}
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
flatten(mesh::PopMesh) = FlatMesh(mesh)
unflatten(mesh::FlatMesh) = mesh.pop_mesh
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
function simple_sigmoid_fn(x, a, theta)
    return @. (1 / (1 + exp(-a * (x - theta))))
end

function sigmoid_fn(x, a, theta)
    return max.(0, simple_sigmoid_fn(x, a, theta) .- simple_sigmoid_fn(0, a, theta))
end
function make_stimulus_fn(mesh; name=nothing, args...)
    stimulus_factories = Dict(
        "smooth_bump" => smooth_bump_factory,
        "sharp_bump" => sharp_bump_factory
    )
    return stimulus_factories[name](mesh; args...)
end
function make_smooth_bump_frame(mesh_coords::Array{DistType}, width::DistType, strength::NumType, steepness::NumType)
    @. strength * (simple_sigmoid_fn(mesh_coords, steepness, -width/2) - simple_sigmoid_fn(mesh_coords, steepness, width/2))
end

function smooth_bump_factory(mesh::AbstractMesh;
                             width=nothing, strength=nothing, duration=nothing,
                             steepness=nothing)
    # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_smooth_bump_frame(coords(mesh), width, strength, steepness)
    return (t) -> @. on_frame * (1 - simple_sigmoid_fn(t, steepness, duration))
end
function make_sharp_bump_frame(mesh::PopMesh, width::DistType, strength::NumType)
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
function sharp_bump_factory(mesh; width=nothing, strength=nothing, duration=nothing)
        # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_sharp_bump_frame(mesh, width, strength)
    off_frame = zeros(on_frame)
    return (t) -> (t <= duration) ? on_frame : off_frame
end
function distance_matrix(xs::SpaceDim)
    # aka Hankel, but that method isn't working in SpecialMatrices
    distance_mx = zeros(eltype(xs), length(xs), length(xs))
    for i in range(1, length(xs))
        distance_mx[:, i] = abs.(xs - xs[i])
    end
    return distance_mx'
end
function sholl_matrix(amplitude::NumType, spread::NumType,
                      dist_mx::Array{NumType,2}, step_size::NumType)
    conn_mx = @. amplitude * step_size * exp(
        -abs(dist_mx / spread)
    ) / (2 * spread)
    return conn_mx
end
function sholl_connectivity(mesh::PopMesh, W::Array{NumType,2}, Σ::Array{NumType,2})::InteractionTensor
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
function sholl_connectivity(mesh::FlatMesh, args...)
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


WilsonCowan73!(dA,A::SpaceState1D,p::WilsonCowan73Params{InteractionTensor},t) = begin
    # Use dA as intermediate variable for tensor op since it is preallocated
    @tensor dA[x_tgt, pop_tgt] = p.W[x_tgt, pop_tgt, x_src, pop_src] * A[x_src, pop_src]
    dA .= (-p.α .* A + p.β .* (1 .- A) .* sigmoid_fn(dA + p.stimulus_fn(t), p.a, p.θ)) ./ p.τ
end

WilsonCowan73!(dA,A::SpaceState1DFlat,p::WilsonCowan73Params{Interaction1DFlat},t) = begin
    # Use dA as intermediate variable for tensor op since it is preallocated
    dA .= (-p.α .* A + p.β .* (1 .- A) .* sigmoid_fn(p.W*A + p.stimulus_fn(t), p.a, p.θ)) ./ p.τ
end

function solve_WilsonCowan73(; model=nothing, solver=nothing, other...)
    solver_params = deepcopy(Dict{Any,Any}(solver))
    model_params = WilsonCowan73Params(model)

    u0 = zeros(model_params.mesh)
    tspan = (0.0, pop!(solver_params, :T))
    prob::DEProblem = ODEProblem(WilsonCowan73!, u0, tspan, model_params)

    if :dt in keys(solver_params)
        solver_params[:alg] = Euler()
        solver_params[:adaptive] = false
    elseif :stiff in keys(solver_params)
        if pop!(solver_params, :stiff) > 0
            solver_params[:alg_hints] = [:stiff]
        end
    end

    if :stiff in keys(solver_params)
        error("Incompatible solver parameters.")
    end

    soln::DESolution = solve(prob; solver_params...)

    return soln
end


macro safe_write(path, writer)
    quote
	if !(isfile($(esc(path))))
	    $(esc(writer))
	else
	    warn("Tried to write existing file: $(esc(path))")
	end
    end
end

function output_dir_name(; root=nothing, simulation_name=nothing, other...)
    now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
    dir_name = joinpath(root, simulation_name, now)
    mkpath(dir_name)
    return dir_name
end

function write_params(dir_name; params...)
    save_path = joinpath(dir_name, "parameters.json")
    @safe_write(save_path, write(save_path, JSON.json(params)))
end


function standardize_timeseries(timeseries, mesh::AbstractMesh)
    # Join array of arrays into matrix Other Dims x Time
    cat(true_ndims(mesh)+1, [standardize_frame(frame, mesh) for frame in timeseries]...)
end
function standardize_frame(frame, mesh::FlatMesh)
    reshape(frame, size(mesh.pop_mesh))
end
function standardize_frame(frame, mesh::PopMesh)
    frame # The PopMesh shape is the standard.
end


  function solution_gif(t, timeseries::Array{NumType,3}; dir_name="", file_name="solution.gif",
                        disable=0, subsample=1, fps=15)
      if disable != 0
          return
      end
      max_activity = maximum(timeseries, (1,2,3))[1] # I don't know why this index is here.
      min_activity = minimum(timeseries, (1,2,3))[1]
      subsample = floor(Int, subsample)
      anim = @animate for i in 1:subsample:length(t)
          plot([timeseries[:,1,i], timeseries[:,2,i]], ylim=(min_activity, max_activity), title="t=$(t[i])")
      end
      save_path = joinpath(dir_name, file_name)
      @safe_write(save_path, gif(anim, save_path, fps=fps))
  end


  function analyse_WilsonCowan73_solution(soln; analyses=nothing, all_params...)
      dir_name = output_dir_name(; analyses...)
      write_params(dir_name; analyses=analyses, all_params...)
      timeseries = standardize_timeseries(soln.u, soln.prob.p.mesh)
      solution_gif(soln.t, timeseries; dir_name=dir_name, analyses[:activity_gif]...)
  end

  function run_WilsonCowan73_trial(json_filename::String, modifications=nothing::Union{Dict, Void})
      all_params = load_WilsonCowan73_parameters(json_filename, modifications)
      solution = solve_WilsonCowan73(; all_params...)
      analyse_WilsonCowan73_solution(solution; all_params...)
  end

export run_WilsonCowan73_trial
end
