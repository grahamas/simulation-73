"A AbstractModel specifies all parameters of a system."
abstract type AbstractModel{T,N,P} <: AbstractParameter{T} end
abstract type AbstractModelwithDelay{T,N,P} <: AbstractModel{T,N,P} end

# FIXME not real dispatch, since it's just an alias
@inline population(A::AbstractArray{T,N}, i) where {T,N} = view_slice_last(A, i)
function population_coordinates(coordinates::AbstractArray{<:CartesianIndex,N}, P) where N
    cat([[CartesianIndex(coord,i) for coord in coordinates] for i in 1:P]...; dims=N+1)
end
population_repeat(arr::AbstractArray{T,N}, P) where {T,N} = repeat(arr, outer=([1 for _ in 1:N]..., P))
"""
    population_timepoint(solution, pop_dx, time_dx)

Return spatial frame for a given population `pop_dx` and time `time_dx`.
"""
@generated function population_timepoint(solution::DiffEqBase.AbstractODESolution{T,NPT}, pop_dx::Int, time_dx::Int) where {T,NPT}
    N = NPT - 2 # N + pop + time
    colons = [:(:) for i in 1:N]
    :(solution[$(colons...), pop_dx, time_dx])
end
# @inline population(A::AbstractArray{T,N}, i) where {T,N} = view_slice_first(A, i)
# function population_coordinates(coordinates::AbstractArray{CartesianIndex,N}, P) where N
#     cat(reshape([CartesianIndex(i,coord) for coord in coordinates], (1, size(coordinates...))) for i in 1:P, dims=1)
# end
# population_repeat(arr::AbstractArray{T,N}, P) where {T,N} = repeat(reshape(arr, (1, size(arr)...)), outer=(P, [1 for _ in 1:N]...))
# """
#     population_timepoint(solution, pop_dx, time_dx)
#
# Return spatial frame for a given population `pop_dx` and time `time_dx`.
# """
# @generated function population_timepoint(solution::DiffEqBase.AbstractODESolution{T,NPT}, pop_dx::Int, time_dx::Int) where {T,NPT}
#     N = NPT - 2 # N + pop + time
#     colons = [:(:) for i in 1:N]
#     :(solution[pop_dx, $(colons...), time_dx])
# end

initial_value(::AbstractModel{T,N,P}, space::AbstractSpace{T,N}) where {T,N,P} = population_repeat(zeros(T,size(space)...), P)

"A Simulation holds an AbstractModel to be solved, the space on which to solve it, the time for which to solve it, the initial value, and various solver options."
struct Simulation{T,M<:AbstractModel{T},S<:AbstractSpace{T}} <: AbstractParameter{T}
    model::M
    space::S
    tspan::Tuple{T,T}
    initial_value::AbstractArray{T}
    solver_options
end
function Simulation(model::M; space::S, tspan, initial_value=initial_value(model, space), opts...) where {T,N,P,M<:AbstractModel{T,N,P}, S<:AbstractSpace{T,N}}
    Simulation{T,M,S}(model, space, tspan, initial_value, opts)
end
# function DrWatson.default_prefix(s::Simulation)
#     space_dir = default_prefix(s.space)
#     nonl_dir = default_prefix(s.nonlinearity)
#     conn_dir = default_prefix(s.connectivity)
#     stim_dir = default_prefix(s.stimulus)
#     joinpath(space_dir, nonl_dir, conn_dir, stim_dir)
# end
# combine_names(name::Symbol, subname::Symbol) = (name, subname)
# combine_names(name::Symbol, subnames::Tuple) = (name, subnames...)
# function DrWatson.allaccess(ap::AbstractParameter)
#     names = fieldnames(typeof(ap))
#     all_names = []
#     for name ∈ names
#         @show typeof(getproperty(ap, name))
#         subnames =  DrWatson.allaccess(getproperty(ap, name))
#         if length(subnames) == 0
#             push!(all_names, name)
#         else
#             push!(all_names, combine_names.(name, subnames)...)
#         end
#     end
#     return all_names
# end


"An Execution holds a Simulation and the solution obtained by running the Simulation."
struct Execution{T,S<:Simulation{T},D<:DESolution}
    simulation::S
    solution::D
end
Execution(s::S) where {T,S <: Simulation{T}} = Execution(s,solve(s))
# DrWatson.savename(e::Execution) = savename(e.simulation)

execute(s::Simulation) = Execution(s)

# initial_value(sim::Simulation) = initial_value(sim.solver)
# time_span(sim::Simulation) = time_span(sim.solver)
# history(simulation::Simulation) = history(simulation.solver)
coordinates(sim::Simulation) = coordinates(sim.space)
coordinates(ex::Execution) = coordinates(ex.simulation)
timepoints(ex::Execution) = ex.solution.t
origin_idx(sim::Simulation) = origin_idx(sim.space)
origin_idx(ex::Execution) = origin_idx(ex.simulation)
# saved_dt(sim::Simulation{T}) where T = saved_dt(sim.solver)
# saved_dx(sim::Simulation{T}) where T = saved_dx(sim.model, sim.solver)
# space_index_info(sim::Simulation{T}) where T = space_index_info(sim.solver)
# time_index_info(sim::Simulation{T}) where T = get_time_index_info(sim.solver)
#
#
#
# function subsampling_idxs(simulation::Simulation{T,<:AbstractModel{T}}, space_subsampler::Subsampler, time_subsampler::Subsampler) where T
#     subsampling_idxs(simulation.model, simulation.solver, space_subsampler, time_subsampler)
# end
#
# function subsample(execution::Execution{T,<:Simulation{T,<:AbstractModel{T}}}; time_subsampler, space_subsampler) where T
#     simulation = execution.simulation
#     solution = execution.solution
#
#     t = timepoints(simulation)
#     x = coordinates(simulation)
#
#     x_dxs, pop_dxs, t_dxs = subsampling_idxs(simulation, space_subsampler, time_subsampler)
#
#     t = t[t_dxs]
#     x = x[x_dxs] # TODO: Remove 1D return assumption
#     wave = solution[x_dxs,pop_dxs,t_dxs]
#
#     return (t,x,wave)
# end


"""
    make_system_mutator(simulation)

Construct the differential function to be provided to the ODE solver.
"""
function make_system_mutator(sim::SIM) where {SIM <: Simulation}
    make_system_mutator(sim.model, sim.space)
end

"""
    generate_problem(simulation)

Return an ODEProblem of the `simulation.model` with time span specified by `simulation.solver`.
"""
function generate_problem(simulation::Simulation{T}) where {T}
    system_mutator! = make_system_mutator(simulation)
    ode_fn = convert(ODEFunction{true}, system_mutator!)
    return ODEProblem(ode_fn, simulation.initial_value, simulation.tspan)
end

# TODO: Add history functionality
# function generate_problem(simulation::Simulation{T,MwD}) where {T, MwD<:AbstractModelwithDelay}
#     system_mutator! = make_system_mutator(simulation)
#     return DDEProblem(system_mutator!, simulation.initial_value, history(simulation), simulation.tspan)
# end
parse_save_idxs(simulation::Simulation, save_idx::Union{Number,AbstractArray}) = save_idx
function parse_save_idxs(simulation::Simulation{T,M}, subsampler::AbstractSubsampler) where {T,N,P,M<:AbstractModel{T,N,P}}
	one_pop_coordinates = coordinate_indices(simulation.space, subsampler)
	population_coordinates(one_pop_coordinates, P)
end

function unpacking_solve(simulation::Simulation; save_idxs=nothing, solver_options...)
    if save_idxs != nothing
        solver_options = (solver_options..., save_idxs=parse_save_idxs(simulation,save_idxs))
    end
    problem = generate_problem(simulation)
    solve(problem; solver_options...)
end

function solve(simulation::Simulation)
    unpacking_solve(simulation; simulation.solver_options...)
end

"""
    run_simulation(jl_filename)

Loads a simulation object defined in `jl_filename`, and save the parameters.
"""
function run_simulation(jl_filename::AbstractString)
    include(jl_filename)
    filecopy(simulation.output, jl_filename, basename(jl_filename))
    @assert all([@isdefined(simulation), @isdefined(output), @isdefined(analysis)])
    execution = execute(simulation)
    results = analyse.(analyses, Ref(execution), Ref(output))
    return (execution, results)
end
