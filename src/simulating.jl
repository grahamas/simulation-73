"A AbstractModel specifies all parameters of a system."
abstract type AbstractModel{T,N,P} <: AbstractParameter{T} end
abstract type AbstractModelwithDelay{T,N,P} <: AbstractModel{T,N,P} end

n_populations(::AbstractModel{T,N,P}) where {T,N,P} = P

initial_value(::AbstractModel{T,N,P}, space::AbstractSpace{T,N}) where {T,N,P} = population_repeat(zeros(space), P)

"A Simulation holds an AbstractModel to be solved, the space on which to solve it, the time for which to solve it, the initial value, and various solver options."
struct Simulation{T,M<:AbstractModel{T},S<:AbstractSpace{T}} <: AbstractParameter{T}
    model::M
    space::S
    tspan::Tuple{T,T}
    initial_value::AbstractArray{T}
    algorithm
    dt::Union{T,Nothing}
    solver_options
end
function Simulation(model::M; space::S, tspan, initial_value=initial_value(model,space), dt=nothing, algorithm, opts...) where {T,N,P,M<:AbstractModel{T,N,P}, S<:AbstractSpace{T,N}}
    Simulation{T,M,S}(model, space, tspan, initial_value, algorithm, dt, opts)
end


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
_space(sp::AbstractSpace, opts) = subsample(sp, get(opts, :save_idxs, nothing))
space(sim::Simulation) = _space(sim.space, sim.solver_options)
space(ex::Execution) = space(ex.simulation)
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
    system_fn! = simulation.model(simulation.space)
    ode_fn = convert(ODEFunction{true}, system_fn!)
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

function unpacking_solve(simulation::Simulation, alg; save_idxs=nothing, dt, solver_options...)
    if save_idxs != nothing
        solver_options = (solver_options..., save_idxs=parse_save_idxs(simulation,save_idxs))
    end
    if dt != nothing
        solver_options = (solver_options..., dt=dt)
    end
    problem = generate_problem(simulation)
    solve(problem, alg; solver_options...)
end

function solve(simulation::Simulation)
    unpacking_solve(simulation, simulation.algorithm; dt=simulation.dt, simulation.solver_options...)
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
