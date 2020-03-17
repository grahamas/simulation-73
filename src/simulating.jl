"A AbstractModel specifies all parameters of a system."
abstract type AbstractModel{T,N,P} <: AbstractParameter{T} end
abstract type AbstractODEModel{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractModelwithDelay{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractNoisyModel{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractNoisySpaceAction{T,N} <: AbstractSpaceAction{T,N} end
abstract type AbstractSimulation{T} <: AbstractParameter{T} end 
abstract type AbstractExecution{T} end 
abstract type AbstractFullExecution{T,SIM} <: AbstractExecution{T} end
export AbstractFullExecution


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

export AbstractNoisyModel, NoisyInputModel, NoisyInputAction, AbstractODEModel

struct NoisyInputModel{T,N,P,M<:AbstractModel{T,N,P}} <: AbstractNoisyModel{T,N,P}
    model::M
    noise_process::NoiseProcess
end
function Base.getproperty(wnm::NoisyInputModel, sym::Symbol)
    if sym ∈ [:noise_process, :model]
        return getfield(wnm, sym)
    else
        return getproperty(getfield(wnm, :model), sym)
    end
end

struct NoisyInputAction{T,N,SA<:AbstractSpaceAction{T,N}} <: AbstractNoisySpaceAction{T,N}
    space_action::SA
    noise_process::NoiseProcess
end
function (act::NoisyInputAction)(du, u, p, t, W)
    du .= W
    act.space_action(du, u, p, t)
end


function (nim::NoisyInputModel)(args...) 
    inner_fn = nim.model(args...)
    NoisyInputAction(inner_fn, nim.noise_process)
end
    

n_populations(::AbstractModel{T,N,P}) where {T,N,P} = P

initial_value(::AbstractModel{T,N,P}, space::AbstractSpace{T,N}) where {T,N,P} = population_repeat(zeros(space), P)

"A Simulation holds an AbstractModel to be solved, the space on which to solve it, the time for which to solve it, the initial value, and various solver options."
struct Simulation{T,M<:AbstractModel{T},S<:AbstractSpace{T},SV<:Union{Tuple{Function,DataType},Nothing}} <: AbstractSimulation{T}
    model::M
    space::S
    tspan::Tuple{T,T}
    initial_value::AbstractArray{T}
    algorithm
    dt::Union{T,Nothing}
    step_reduction::SV
    global_reduction::Function
    solver_options
end
function Simulation(model::M; space::S, tspan, initial_value=initial_value(model,space), algorithm, dt=nothing, step_reduction::SV=nothing, global_reduction=identity, opts...) where {T,N,P,M<:AbstractModel{T,N,P}, S<:AbstractSpace{T,N},SV}
    return Simulation{T,M,S,SV}(model, space, tspan, initial_value, algorithm, dt, step_reduction, global_reduction, opts)
end
function Simulation(model::Missing; kwargs...)
    return FailedSimulation{Missing}()
end

"An Execution holds a Simulation and the solution obtained by running the Simulation."
struct Execution{T,S<:AbstractSimulation{T},D<:DESolution} <: AbstractFullExecution{T,S}
    simulation::S
    solution::D
end
struct ReducedExecution{T,S<:AbstractSimulation{T},SV<:SavedValues} <: AbstractExecution{T}
    simulation::S
    saved_values::SV
end
struct AugmentedExecution{T,S<:AbstractSimulation{T},D<:DESolution,SV<:SavedValues} <: AbstractFullExecution{T,S}
    simulation::S
    solution::D
    saved_values::SV
end
make_execution(s::S, sol::DESolution) where {T,S <: Simulation{T}} = Execution(s,sol)
make_execution(s::S, sv::SavedValues) where {T,S <: Simulation{T}} = ReducedExecution(s,sv)
make_execution(s::S, (sol,sv)::Tuple{<:DESolution,<:SavedValues}) where {T,S <: Simulation{T}} = AugmentedExecution(s,sol,sv)
export ReducedExecution, AugmentedExecution,make_execution
function execute(s::Simulation)
    return make_execution(s, solve(s))
end
function execute(s::FailedSimulation)
    return FailedExecution(s)
end
            
coordinates(sim::Simulation) = coordinates(space(sim))
coordinates(ex::AbstractExecution) = coordinates(space(ex))
timepoints(ex::AbstractFullExecution) = ex.solution.t
_space(sp::AbstractSpace, opts) = subsample(sp, get(opts, :save_idxs, nothing))
space(sim::Simulation) = _space(sim.space, sim.solver_options)
space(ex::AbstractExecution) = space(ex.simulation)
origin_idx(sim::Simulation) = origin_idx(sim.space)
origin_idx(ex::AbstractExecution) = origin_idx(ex.simulation)
extrema(ex::AbstractFullExecution) = extrema(ex.solution.u)
stimulus_center(mod::AbstractModel) = center(mod.stimulus)
stimulus_center(sim::Simulation) = stimulus_center(sim.model)
export stimulus_center

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
function generate_problem(simulation::Simulation{T,<:AbstractODEModel}; callback=nothing) where {T}
    system_fn! = simulation.model(simulation.space)
    ode_fn = convert(ODEFunction{true}, system_fn!)
    return ODEProblem(ode_fn, simulation.initial_value, simulation.tspan, callback=callback)
end

function generate_problem(simulation::Simulation{T,<:AbstractNoisyModel}; callback=nothing) where {T}
    system_fn! = simulation.model(simulation.space)
    return RODEProblem(system_fn!, simulation.initial_value, simulation.tspan, noise=simulation.model.noise_process, noise_prototype=zeros(size(simulation.initial_value)...), callback=callback)
end

# TODO: Add history functionality
# function generate_problem(simulation::Simulation{T,MwD}) where {T, MwD<:AbstractModelwithDelay}
#     system_mutator! = make_system_mutator(simulation)
#     return DDEProblem(system_mutator!, simulation.initial_value, history(simulation), simulation.tspan)
# end
parse_save_idxs(simulation::Simulation, save_idx::Union{Number,AbstractArray{<:Number}}) = save_idx
function parse_save_idxs(simulation::Simulation{T,M}, subsampler::Union{AbstractSubsampler,AbstractArray{<:AbstractSubsampler}}) where {T,N,P,M<:AbstractModel{T,N,P}}
	one_pop_coordinates = coordinate_indices(simulation.space, subsampler)
	population_coordinates(one_pop_coordinates, P)
end

function solve(simulation::Simulation, alg; callback=nothing, save_idxs=nothing, dt=nothing, solver_options...)
    if save_idxs != nothing
        solver_options = (solver_options..., save_idxs=parse_save_idxs(simulation, save_idxs))
    end
    if dt != nothing
        solver_options = (solver_options..., dt=dt)
    end
    problem = generate_problem(simulation; callback=callback)
    solve(problem, alg; solver_options...)
end

function solve(simulation::Simulation{A,B,C,Nothing}; callback=nothing, solver_opts...) where {A,B,C}
    if :callback ∈ keys(simulation.solver_options)
        callback = CallbackSet(callback, simulation.solver_options[:callback])
    end
    
    sol = solve(simulation, simulation.algorithm; dt=simulation.dt, simulation.solver_options..., solver_opts..., callback=callback)
    return sol
end

function solve(simulation::Simulation{T,A,B,<:Tuple}; callback=nothing, solver_opts...) where {T,A,B}
    save_func, save_type = simulation.step_reduction
    sv = SavedValues(T,save_type)
    all_callbacks = SavingCallback(save_func, sv)
        
    if :callback ∈ keys(simulation.solver_options)
        all_callbacks = CallbackSet(simulation.solver_options[:callback], all_callbacks)
    end
    
    if callback !== nothing
        all_callbacks = CallbackSet(callback, all_callbacks)
    end
    # save_callback has all callbacks
    
    sol = solve(simulation, simulation.algorithm; dt=simulation.dt, simulation.solver_options..., solver_opts..., callback=all_callbacks)
    return (sol, sv)
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

