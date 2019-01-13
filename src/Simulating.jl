module Simulating

using Modeling
using Records
using Parameters
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve, Euler, ODEProblem
using Meshes

abstract type AbstractPlotSpecification end
@with_kw struct Analyses{T}
    plot_specs::Array{AbstractPlotSpecification}
end
export AbstractPlotSpecification, Analyses

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
save_dt(s::Solver{T}) where T = s.simulated_dt / s.time_save_every

function save_idxs(solver::Solver{T}, space::SP) where {T,P, SP <: PopSpace{T,1,P}}#::Union{Nothing,Array{CartesianIndex}} 
    if solver.space_save_every == 1
        return nothing
    end
    return CartesianIndex.(Iterators.product(axes(space)...))[1:solver.space_save_every:end,:]
end

"""A Simulation object runs its own simulation upon initialization."""
struct Simulation{T,M<:Model{T},S<:Solver{T}}
    model::M
    solver::S
    analyses::Analyses{T}
    output::AbstractOutput
    solution::DESolution
    Simulation{T,M,S}(m,s,a,o) where {T,M<:Model{T},S<:Solver{T}} = new(m,s,a,o,_solve(m,s))
end

function Simulation(; model::M, solver::S, analyses::Analyses{T}, output::AbstractOutput) where {T, M<:Model{T}, S<:Solver{T}}
    Simulation{T,M,S}(model,solver,analyses,output)
end

Modeling.initial_value(sim::Simulation) = initial_value(sim.model)
time_span(solver::Solver{T}) where T = solver.tspan
time_span(sim::Simulation) = time_span(sim.solver)
Modeling.space_arr(sim::Simulation) = space_arr(sim.model)[1:sim.solver.space_save_every:end,1] # TODO: remove 1D assumption
time_arr(sim::Simulation) = sim.solution.t
function Meshes.get_origin(sim::Simulation) # TODO: Remove 1D assumption
    round(Int, Meshes.get_origin(sim.model.space)[1] / sim.solver.space_save_every) 
end

Base.minimum(sim::Simulation) = minimum(map(minimum, sim.solution.u))
Base.maximum(sim::Simulation) = maximum(map(maximum, sim.solution.u))

function write_params(sim::Simulation)
    write_object(sim.output, "parameters.jld2", "sim", sim)
end

"""
    _solve wraps the DifferentialEquations function, solve.
    Note that the method accepting a Simulation object should take a
    partially initialized Simulation.
"""
function _solve(model,solver)
    problem = generate_problem(model, solver)
    _solve(problem, solver, model.space)
end
function _solve(problem::ODEProblem, solver::Solver{T,Euler}, space::PopSpace{T}) where T
    # TODO: Calculate save_idxs ACCOUNTING FOR pops
    solve(problem, Euler(), dt=solver.simulated_dt,
            timeseries_steps=solver.time_save_every,
            save_idxs=save_idxs(solver, space))
end
function _solve(problem::ODEProblem, solver::Solver{T,Nothing}, space::PopSpace{T}) where T
    solve(problem, saveat=save_dt(solver), timeseries_steps=solver.time_save_every,
        save_idxs=save_idxs(solver, space), alg_hints=[solver.stiffness])
end

""" run_simulation loads a simulation object defined in a jl script, and saves the parameters. """
function run_simulation(jl_filename::AbstractString)
    include(jl_filename)
    filecopy(simulation.output, jl_filename, basename(jl_filename))
    return simulation
end

export run_simulation, Simulation, write_params, Solver,
    time_span, time_arr

end