"A AbstractModel specifies all parameters of a system."
abstract type AbstractModel{T,N,P} <: AbstractParameter{T} end
abstract type AbstractODEModel{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractModelwithDelay{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractNoisyModel{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractSimulation{T} <: AbstractParameter{T} end 
abstract type AbstractExecution{T} end 

function (am::Type{<:AbstractModel})(fallback_args...; fallback_kwargs...)
    #@warn """Model $am undefined!
    #---------------------
    #$fallback_args
   # 
    #$fallback_kwargs
    #---------------------   
    #"""
    missing
end

struct FailedSimulation{T} <: AbstractSimulation{T} end
struct FailedExecution{T,S<:AbstractSimulation{T}} <: AbstractExecution{T}
    sim::S
end

export AbstractNoisyModel, WeinerNoisyModel, AbstractODEModel

struct WeinerNoisyModel{T,N,P,M<:AbstractModel{T,N,P}} <: AbstractNoisyModel{T,N,P}
    noise_coefficient::Float64
    model::M
end
function Base.getproperty(wnm::WeinerNoisyModel, sym::Symbol)
    if sym ∈ [:noise_coefficient, :model]
        return getfield(wnm, sym)
    else
        return getproperty(getfield(wnm, :model), sym)
    end
end

(wnm::WeinerNoisyModel)(args...) = (wnm.model(args...), (du,u,p,t) -> (du .= wnm.noise_coefficient))

n_populations(::AbstractModel{T,N,P}) where {T,N,P} = P

initial_value(::AbstractModel{T,N,P}, space::AbstractSpace{T,N}) where {T,N,P} = population_repeat(zeros(space), P)

"A Simulation holds an AbstractModel to be solved, the space on which to solve it, the time for which to solve it, the initial value, and various solver options."
struct Simulation{T,M<:AbstractModel{T},S<:AbstractSpace{T}} <: AbstractSimulation{T}
    model::M
    space::S
    tspan::Tuple{T,T}
    initial_value::AbstractArray{T}
    algorithm
    dt::Union{T,Nothing}
    solver_options
end
function Simulation(model::M; space::S, tspan, initial_value=initial_value(model,space), dt=nothing, algorithm, opts...) where {T,N,P,M<:AbstractModel{T,N,P}, S<:AbstractSpace{T,N}}
    return Simulation{T,M,S}(model, space, tspan, initial_value, algorithm, dt, opts)
end
function Simulation(model::Missing; kwargs...)
    return FailedSimulation{Missing}()
end

"An Execution holds a Simulation and the solution obtained by running the Simulation."
struct Execution{T,S<:AbstractSimulation{T},D<:DESolution} <: AbstractExecution{T}
    simulation::S
    solution::D
end
Execution(s::S) where {T,S <: Simulation{T}} = Execution(s,solve(s))

function execute(s::Simulation)
    return Execution(s)
end
function execute(s::FailedSimulation)
    return FailedExecution(s)
end
            
coordinates(sim::Simulation) = coordinates(space(sim))
coordinates(ex::AbstractExecution) = coordinates(space(ex))
timepoints(ex::Execution) = ex.solution.t
_space(sp::AbstractSpace, opts) = subsample(sp, get(opts, :save_idxs, nothing))
space(sim::Simulation) = _space(sim.space, sim.solver_options)
space(ex::Execution) = space(ex.simulation)
origin_idx(sim::Simulation) = origin_idx(sim.space)
origin_idx(ex::Execution) = origin_idx(ex.simulation)
extrema(ex::AbstractExecution) = extrema(ex.solution.u)

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
function generate_problem(simulation::Simulation{T,<:AbstractODEModel}) where {T}
    system_fn! = simulation.model(simulation.space)
    ode_fn = convert(ODEFunction{true}, system_fn!)
    return ODEProblem(ode_fn, simulation.initial_value, simulation.tspan)
end

function generate_problem(simulation::Simulation{T,<:AbstractNoisyModel}) where {T}
    system_fn!, noise_fn! = simulation.model(simulation.space)
    return SDEProblem(system_fn!, noise_fn!, simulation.initial_value, simulation.tspan)
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

# """
#     run_simulation(jl_filename)

# Loads a simulation object defined in `jl_filename`, and save the parameters.
# """
# function run_simulation(jl_filename::AbstractString)
#     include(jl_filename)
#     filecopy(simulation.output, jl_filename, basename(jl_filename))
#     @assert all([@isdefined(simulation), @isdefined(output), @isdefined(analysis)])
#     execution = execute(simulation)
#     results = analyse.(analyses, Ref(execution), Ref(output))
#     return (execution, results)
# end

