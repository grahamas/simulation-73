
#region Model
abstract type Model{T,N,P} <: AbstractParameter{T} end

function space_arr(model::Model)::AbstractArray
    Calculated(model.space).value
end

function initial_value(model::Model{T,N,P}) where {T,N,P}
    space_zeros = zeros(model.space)
end
#endregion

#region Solver
struct Solver{T,ALG<:Union{OrdinaryDiffEqAlgorithm,Nothing},DT<:Union{T,Nothing}}
    tspan::Tuple{T,T}
    algorithm::ALG
    simulated_dt::DT
    time_save_every::Int
    space_save_every::Int # TODO: Remove 1D Assumption
    stiffness::Symbol
end

function Solver(; start_time::T=0.0, stop_time::T, dt::DT, time_save_every::Int=1,
                    space_save_every::Int=1,
                    algorithm::ALG=nothing, stiffness=:auto) where {T, ALG, DT<:Union{T,Nothing}}
    Solver{T,ALG,DT}((start_time, stop_time), algorithm, dt, time_save_every, space_save_every, stiffness)
end
save_dt(s::Solver{T}) where T = s.simulated_dt * s.time_save_every

function save_idxs(solver::Solver{T}, space::SP) where {T,P, SP <: PopSpace{T,1,P}}#::Union{Nothing,Array{CartesianIndex}} 
    if solver.space_save_every == 1
        return nothing
    end
    return CartesianIndex.(Iterators.product(axes(space)...))[1:solver.space_save_every:end,:]
end
#endregion

#region Simulation
"""A Simulation object runs its own simulation upon initialization."""
struct Simulation{T,M<:Model{T},S<:Solver{T}}
    model::M
    solver::S
    solution::DESolution
    Simulation{T,M,S}(m,s) where {T,M<:Model{T},S<:Solver{T}} = new(m,s,_solve(m,s))
end

function Simulation(; model::M, solver::S) where {T, M<:Model{T}, S<:Solver{T}}
    Simulation{T,M,S}(model,solver)
end

initial_value(sim::Simulation) = initial_value(sim.model)
time_span(solver::Solver{T}) where T = solver.tspan
time_span(sim::Simulation) = time_span(sim.solver)
space_arr(sim::Simulation) = space_arr(sim.model)[1:sim.solver.space_save_every:end,1] # TODO: remove 1D assumption
time_arr(sim::Simulation) = sim.solution.t
function get_origin(sim::Simulation) # TODO: Remove 1D assumption
    round(Int, get_origin(sim.model.space)[1] / sim.solver.space_save_every) 
end
save_dt(sim::Simulation{T}) where T = save_dt(sim.solver)
save_dx(sim::Simulation{T}) where T = step(sim.model.space) * sim.solver.space_save_every

Base.minimum(sim::Simulation) = minimum(map(minimum, sim.solution.u))
Base.maximum(sim::Simulation) = maximum(map(maximum, sim.solution.u))

function write_params(sim::Simulation)
    write_object(sim.output, "parameters.jld2", "sim", sim)
end
#endregion

"""
    _solve wraps the DifferentialEquations function, solve.
    Note that the method accepting a Simulation object should take a
    partially initialized Simulation.
"""
generate_problem() = error("undefined")
function _solve(model,solver)
    problem = generate_problem(model, solver)
    _solve(problem, solver, model.space)
end
function _solve(problem::ODEProblem, solver::Solver{T,Euler}, space::PopSpace{T}) where T
    # TODO: Calculate save_idxs ACCOUNTING FOR pops
    @show "Solving Euler"
    solve(problem, Euler(), dt=solver.simulated_dt,
            #saveat=save_dt(solver),
            timeseries_steps=solver.time_save_every,
            save_idxs=save_idxs(solver, space))
end
function _solve(problem::ODEProblem, solver::Solver{T,Nothing}, space::PopSpace{T}) where T
    @show "Solving with default ALG"
    solve(problem, saveat=save_dt(solver), timeseries_steps=solver.time_save_every,
        save_idxs=save_idxs(solver, space), alg_hints=[solver.stiffness])
end

""" run_simulation loads a simulation object defined in a jl script, and saves the parameters. """
function run_simulation(jl_filename::AbstractString)
    include(jl_filename)
    filecopy(simulation.output, jl_filename, basename(jl_filename))
    return simulation
end

function generate_problem(model::M, solver::SV) where {T,M<:Model{T},SV<:Solver{T}}
    tspan = time_span(solver)
    u0 = initial_value(model)

    calculated_model = Calculated(model)

    system_fn! = make_calculated_function(calculated_model)

    ode_fn = convert(ODEFunction{true}, system_fn!)
    return ODEProblem(ode_fn, u0, tspan, nothing)
end