module Simulating

using Modeling
using Analysis
using Records
using Parameters
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve, Euler, ODEProblem


abstract type AbstractSolver{T} end


@with_kw struct AutoSolver{S} <: AbstractSolver{S}
    T::S
end
@with_kw struct EulerSolver{S} <: AbstractSolver{S}
    T::S
    dt::S
end

time_span(solver::AbstractSolver{S}) where S = (zero(S), solver.T)

struct Simulation{T,M<:Model{T},S<:AbstractSolver{T}}
    model::M
    solver::S
    analyses::Analyses{T}
    output::AbstractOutput
end

function Simulation(; model::M, solver::S, analyses::Analyses{T}, output::AbstractOutput) where {T, M<:Model{T}, S<:AbstractSolver{T}}
    Simulation{T,M,S}(model,solver,analyses,output)
end

Modeling.initial_value(sim::Simulation) = initial_value(sim.model)
time_span(sim::Simulation) = time_span(sim.solver)
solver_params(sim::Simulation) = sim.solver.params

function write_params(sim::Simulation)
	write_object(sim.output, "parameters.jld2", "sim", sim)
end

function solve(problem::ODEProblem, solver::EulerSolver)
    solve(problem, Euler(), dt=solver.dt)
end
function solve(problem::ODEProblem, solver::AutoSolver)
    solve(problem)
end

function simulate(jl_filename::AbstractString)
	include(jl_filename)
	#@load "parameters.jld2" simulation
    filecopy(simulation.output, jl_filename, basename(jl_filename))
	simulate(simulation)
end

function simulate(simulation::Simulation)
    problem = generate_problem(simulation)
	solution = solve(problem, simulation.solver)
	analyse(simulation, solution)
end

function Analysis.analyse(simulation::Simulation, solution::DESolution)
    analyses = simulation.analyses
    results = Results(simulation.model, solution, analyses.subsampler)
    analyse(analyses, results, simulation.output)
end

function results_only(jl_filename::AbstractString)
    include(jl_filename)
    problem = generate_problem(simulation)
    solution = solve(problem, simulation.solver)
    return (simulation, Results(simulation.model, solution, simulation.analyses.subsampler))
end

export simulate, Simulation, write_params, Solver,
    time_span, solve, AbstractSolver

end