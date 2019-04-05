"A Model specifies all parameters of a system."
abstract type Model{T,N,P} <: AbstractParameter{T} end

"A Solver holds the parameters of the differential equation solver."
struct Solver{T,ALG<:Union{OrdinaryDiffEqAlgorithm,Nothing},DT<:Union{T,Nothing}}
    tspan::Tuple{T,T}
    algorithm::ALG
    simulated_dt::DT
    time_save_every::Int
    space_save_every::Int # TODO: Remove 1D Assumption
    stiffness::Symbol
    #dense::Bool
end

function Solver{S}(; start_time::S=0.0, stop_time::S, dt::DT, time_save_every::Int=1, space_save_every=1, algorithm::ALG=nothing, stiffness=:auto) where {S, ALG, DT<:Union{S,Nothing}}
    Solver{S,ALG,DT}((start_time, stop_time), algorithm, dt, time_save_every, space_save_every, stiffness)
end
"""
    save_dt(solver)

Return the time step saved by a solver.
"""
save_dt(s::Solver{T}) where T = s.simulated_dt * s.time_save_every

"""
    save_idxs(solver, space)

Return the indices into space of the values that the solver saves.
"""
function save_idxs(solver::Solver{T}, space::SP) where {T,P, SP <: Pops{P,T}}#::Unio{P,n}{Nothing,Array{CartesianIndex}}
    if all(solver.space_save_every .== 1)
        return nothing
    end
    all_indices = CartesianIndices(space)
    space_saved_subsample(all_indices, solver)
end

"A Simulation holds a Model to be solved and a Solver to solve it."
struct Simulation{T,M<:Model{T},S<:Solver{T}}
    model::M
    solver::S
    Simulation{T,M,S}(m::M,s::S) where {T,M<:Model{T},S<:Solver{T}} = new(m,s)
end
function Simulation(; model::M, solver::S) where {T, M<:Model{T}, S<:Solver{T}}
    Simulation{T,M,S}(model,solver)
end

"An Execution holds a Simulation and the solution obtained by running the Simulation."
struct Execution{T,S<:Simulation{T},D<:DESolution}
    simulation::S
    solution::D
end
Execution(s::S) where {T,S <: Simulation{T}} = Execution(s,solve(s))

"""
    initial_value(model)
    initial_value(simulation)

Return the model's initial value (defaults to all zeros)
"""
initial_value(model::Model{T,N,P}) where {T,N,P} = zero(model.space)
initial_value(sim::Simulation) = initial_value(sim.model)
"""
    time_span(solver)
    time_span(simulation)

Return the time span over which the solver runs.
"""
time_span(solver::Solver{T}) where T = solver.tspan
time_span(sim::Simulation) = time_span(sim.solver)
"""
    space_saved_subsample_(arr, solver)

Return the part of `arr` that would have been saved by `solver`.
"""
function _space_saved_subsample_(arr, solver::Solver)
    collect(arr)[[StrideToEnd(i) for i in solver.space_save_every]...,:]
end
"""
    saved_space_arr(model, solver)
    saved_space_arr(simulation)

Return the spatial coordinates of values saved by `solver`
"""
saved_space_arr(model::Model, solver::Solver) = space_saved_subsample_(coordinates(model.space), solver)
saved_space_arr(sim::Simulation) = saved_space_arr(sim.model, sim.solver)
"""
    saved_time_arr(solver)
    saved_time_arr(simulation)

Return the times saved by `solver`.
"""
function saved_time_arr(solver::Solver{T}) where T
    start, stop = time_span(solver)
    start:save_dt(solver):stop
end
saved_time_arr(sim::Simulation) = saved_time_arr(sim.solver)#sim.solution.t

"""
    get_space_origin_idx(model)
    get_space_origin_idx(model, solver)
    get_space_origin_idx(simulation)

Return the index of the spatial origin of `model`'s `space`.
"""
function get_space_origin_idx(model::Model)
    get_space_origin_idx(model.space)
end
function get_space_origin_idx(model::Model, solver::Solver)  # TODO: Remove 1D assumption
    round(Int, get_space_origin_idx(model) ./ solver.space_save_every)
end
function get_space_origin_idx(sim::Simulation)
    get_space_origin_idx(sim.model, sim.solver)
end

save_dt(sim::Simulation{T}) where T = save_dt(sim.solver)
save_dx(model::Model, solver::Solver)= step(model.space) * solver.space_save_every
save_dx(sim::Simulation{T}) where T = save_dx(sim.model, sim.solver)
Base.minimum(solution::DESolution) = minimum(map(minimum, solution.u))
Base.maximum(solution::DESolution) = maximum(map(maximum, solution.u))

"""
    get_space_index_info(model, solver)
    get_space_index_info(simulation)

Return IndexInfo for saved space array.
"""
get_space_index_info(model::Model{T}, solver::Solver{T}) where T = IndexInfo(save_dx(model, solver), get_space_origin_idx(model, solver))
get_space_index_info(sim::Simulation{T}) where T = get_space_index_info(sim.model, sim.solver)
"""
    get_time_index_info(solver)
    get_time_index_info(simulation)

Return IndexInfo for saved time array.
"""
get_time_index_info(solver::Solver{T}) where T = IndexInfo(save_dt(solver), 1)
get_time_index_info(sim::Simulation{T}) where T = get_time_index_info(sim.solver)

"""
    pop_frame(solution, pop_dx, time_dx)

Return spatial frame for a given population `pop_dx` and time `time_dx`.
"""
@generated function pop_frame(solution::ODESolution{T,NPT,<:Array{<:Array{T,NP},1}}, pop_dx::Int, time_dx::Int) where {T,NP,NPT}
    N = NP - 1
    colons = [:(:) for i in 1:N]
    :(solution[$(colons...),pop_dx, time_dx])
end

"""
    write_params(simulation)

Write the simulation to a top-level parameters file hard-coded "parameters.jld2"
"""
function write_params(sim::Simulation)
    write_object(sim.output, "parameters.jld2", "sim", sim)
end


function subsampling_idxs(model::Model, solver::Solver, time_subsampler, space_subsampler)
    x_info = get_space_index_info(model, solver)
    t_info = get_time_index_info(solver)

    x_dxs = subsampling_idxs(x_info, space_subsampler)
    t_dxs = subsampling_idxs(t_info, time_subsampler)

    return (x_dxs, 1, t_dxs)
end

function subsampling_idxs(simulation::Simulation{T,<:Model{T,1}}, time_subsampler::Subsampler, space_subsampler::Subsampler) where T
    subsampling_idxs(simulation.model, simulation.solver, time_subsampler, space_subsampler)
end
function subsampling_time_idxs(solver::Solver, t_target::AbstractArray)
    t_solver = time_span(solver)[1]:save_dt(solver):time_span(solver)[end]
    subsampling_idxs(t_target, t_solver)
end
function subsampling_space_idxs(model::Model, solver::Solver, x_target::AbstractArray)
    x_model = saved_space_arr(model, solver)
    subsampling_idxs(x_target, x_model)
end

function subsample(execution::Execution{T,<:Simulation{T,<:Model{T}}}; time_subsampler, space_subsampler) where T
    simulation = execution.simulation
    solution = execution.solution

    t = saved_time_arr(simulation)
    x = saved_space_arr(simulation)

    x_dxs, pop_dxs, t_dxs = subsampling_idxs(simulation, time_subsampler, space_subsampler)

    t = t[t_dxs]
    x = x[x_dxs] # TODO: Remove 1D return assumption
    wave = solution[x_dxs,pop_dxs,t_dxs]

    return (t,x,wave)
end


"""
    generate_problem(model, solver)

Return an ODEProblem of the `model` with time span specified by `solver`.
"""
function generate_problem(model::M, solver::SV) where {T,M<:Model{T},SV<:Solver{T}}
    tspan = time_span(solver)
    u0 = initial_value(model)

    calculated_model = Calculated(model)

    system_fn! = make_calculated_function(calculated_model)

    ode_fn = convert(ODEFunction{true}, system_fn!)
    return ODEProblem(ode_fn, u0, tspan, nothing)
end
generate_problem(simulation) = generate_problem(simulation.model, simulation.solver)
solve(simulation::Simulation) = _solve(simulation.model, simulation.solver)
function _solve(model::Model,solver::Solver)
    problem = generate_problem(model, solver)
    _solve(problem, solver, model.space)
end
function _solve(problem::ODEProblem, solver::Solver{T,Euler}, space::Pops{P,T}) where {P,T}
    # TODO: Calculate save_idxs ACCOUNTING FOR pops
    @show "Solving Euler"
    solve(problem, Euler(), dt=solver.simulated_dt,
            #saveat=save_dt(solver),
            timeseries_steps=solver.time_save_every,
            save_idxs=save_idxs(solver, space))
end
function _solve(problem::ODEProblem, solver::Solver{T,Nothing}, space::Pops{P,T}) where {P,T}
    @show "Solving with default ALG"
    solve(problem, saveat=save_dt(solver), timeseries_steps=solver.time_save_every,
        save_idxs=save_idxs(solver, space), alg_hints=[solver.stiffness])
end

"""
    run_simulation(jl_filename)

Loads a simulation object defined in `jl_filename`, and save the parameters. 
"""
function run_simulation(jl_filename::AbstractString)
    include(jl_filename)
    filecopy(simulation.output, jl_filename, basename(jl_filename))
    return execution
end
