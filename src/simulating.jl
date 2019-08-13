"A AbstractModel specifies all parameters of a system."
abstract type AbstractModel{T} <: AbstractParameter{T} end
abstract type AbstractModelwithDelay{T} <: AbstractModel{T} end

"A Simulation holds an AbstractModel to be solved and a Solver to solve it."
struct Simulation{T,M<:AbstractModel{T},S<:Solver{T}} <: AbstractParameter{T}
    model::M
    solver::S
    Simulation{T,M,S}(m::M,s::S) where {T,M<:AbstractModel{T},S<:Solver{T}} = new(m,s)
end
function Simulation(; model::M, solver::S) where {T, M<:AbstractModel{T}, S<:Solver{T}}
    Simulation{T,M,S}(model,solver)
end
function DrWatson.default_prefix(s::Simulation)
    space_dir = default_prefix(s.space)
    nonl_dir = default_prefix(s.nonlinearity)
    conn_dir = default_prefix(s.connectivity)
    stim_dir = default_prefix(s.stimulus)
    joinpath(space_dir, nonl_dir, conn_dir, stim_dir)
end
combine_names(name::Symbol, subname::Symbol) = (name, subname)
combine_names(name::Symbol, subnames::Tuple) = (name, subnames...)
function DrWatson.allaccess(ap::AbstractParameter)
    names = fieldnames(typeof(ap))
    all_names = []
    for name ∈ names
        @show typeof(getproperty(ap, name))
        subnames =  DrWatson.allaccess(getproperty(ap, name))
        if length(subnames) == 0
            push!(all_names, name)
        else
            push!(all_names, combine_names.(name, subnames)...)
        end
    end
    return all_names
end


"An Execution holds a Simulation and the solution obtained by running the Simulation."
struct Execution{T,S<:Simulation{T},D<:DESolution}
    simulation::S
    solution::D
end
Execution(s::S) where {T,S <: Simulation{T}} = Execution(s,solve(s))
DrWatson.savename(e::Execution) = savename(e.simulation)

execute(s::Simulation) = Execution(s)

initial_value(sim::Simulation) = initial_value(sim.solver)
time_span(sim::Simulation) = time_span(sim.solver)
history(simulation::Simulation) = history(simulation.solver)
saved_coordinates(sim::Simulation) = saved_coordinates(sim.solver)
saved_coordinates(ex::Execution) = saved_coordinates(ex.simulation)
timepoints(sim::Simulation) = timepoints(sim.solver)#sim.solution.t
timepoints(ex::Execution) = ex.solution.t
function origin_idx(sim::Simulation)
    origin_idx(sim.solver)
end
origin_idx(ex::Execution) = origin_idx(ex.simulation)
saved_dt(sim::Simulation{T}) where T = saved_dt(sim.solver)
saved_dx(sim::Simulation{T}) where T = saved_dx(sim.model, sim.solver)
space_index_info(sim::Simulation{T}) where T = space_index_info(sim.solver)
time_index_info(sim::Simulation{T}) where T = get_time_index_info(sim.solver)



function subsampling_idxs(simulation::Simulation{T,<:AbstractModel{T}}, space_subsampler::Subsampler, time_subsampler::Subsampler) where T
    subsampling_idxs(simulation.model, simulation.solver, space_subsampler, time_subsampler)
end

function subsample(execution::Execution{T,<:Simulation{T,<:AbstractModel{T}}}; time_subsampler, space_subsampler) where T
    simulation = execution.simulation
    solution = execution.solution

    t = timepoints(simulation)
    x = coordinates(simulation)

    x_dxs, pop_dxs, t_dxs = subsampling_idxs(simulation, space_subsampler, time_subsampler)

    t = t[t_dxs]
    x = x[x_dxs] # TODO: Remove 1D return assumption
    wave = solution[x_dxs,pop_dxs,t_dxs]

    return (t,x,wave)
end


"""
    generate_problem(model, solver)

Return an ODEProblem of the `model` with time span specified by `solver`.
"""
function generate_problem(simulation::Simulation{T}) where {T}
    system_mutator! = make_system_mutator(simulation)
    ode_fn = convert(ODEFunction{true}, system_mutator!)
    return ODEProblem(ode_fn, initial_value(simulation), time_span(simulation))
end

function generate_problem(simulation::Simulation{T,MwD}) where {T, MwD<:AbstractModelwithDelay}
    system_mutator! = make_system_mutator(simulation)
    return DDEProblem(system_mutator!, initial_value(simulation), history(simulation), time_span(simulation))
end

function solve(simulation::Simulation)
    problem = generate_problem(simulation)
    _solve(problem, simulation.solver, simulation.model)
end
function _solve(problem::ODEProblem, solver::Solver{T,Euler}, model::AbstractModel{T}) where {T}
    # TODO: Calculate save_idxs ACCOUNTING FOR pops
    @show "Solving Euler"
    solve(problem, Euler(), dt=solver.simulated_dt,
            #saveat=saved_dt(solver),
            timeseries_steps=solver.time_save_every,
            save_idxs=save_idxs(model, solver))
end
function _solve(problem::ODEProblem, solver::Solver{T,Nothing}, model::AbstractModel{T}) where {T}
    @show "Solving with default ALG"
    solve(problem, saveat=saved_dt(solver), timeseries_steps=solver.time_save_every,
        save_idxs=save_idxs(model, solver), alg_hints=[solver.stiffness])
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

"""
    make_mutators(simulation)

Construct a list functions to successively mutate the system's difference (`dA`).
"""
make_mutators(sim::Simulation) = make_mutators(sim.model)

"""
    make_system_mutator(simulation)

Construct the differential function to be provided to the ODE solver.
"""
function make_system_mutator(sim::SIM) where {SIM <: Simulation}
    make_system_mutator(sim.model)
end
