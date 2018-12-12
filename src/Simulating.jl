module Simulating

using Modeling
using Analysis
using Records
using Parameters
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve, Euler


@with_kw struct Solver{S,ODEA<:Union{OrdinaryDiffEqAlgorithm,Nothing}}
    T::S
    kwargs::Dict
    solution_method::ODEA=nothing
end
time_span(solver::Solver{T}) where T = (zero(T), solver.T)

@with_kw struct Simulation{T,M<:Model{T}}
    model::M
    solver::Solver{T}
    analyses::Analyses{T}
    output::AbstractOutput
end

Modeling.initial_value(sim::Simulation) = initial_value(sim.model)
time_span(sim::Simulation) = time_span(sim.solver)
solver_params(sim::Simulation) = sim.solver.params

function solver_args(solver::Solver{S,Nothing}) where S
    return ((), solver.kwargs)
end

function solver_args(solver::Solver{S,Euler}) where S
    return ((solver.solution_method,), solver.kwargs)
end

function write_params(sim::Simulation)
	write_object(sim.output, "parameters.jld2", "sim", sim)
end

function solve(simulation::Simulation)
    problem = generate_problem(simulation)
    args, kwargs = solver_args(simulation.solver)
    soln = solve(problem, args...; kwargs...)
    return soln
end

function simulate(jl_filename::AbstractString)
	include(jl_filename)
	@load "parameters.jld2" simulation
    filecopy(simulation.output, jl_filename, basename(jl_filename))
	simulate(simulation)
end

function simulate(simulation::Simulation)
	solution = solve(simulation)
	analyse(simulation, solution)
end

function Analysis.analyse(simulation::Simulation, solution::DESolution)
    analyses = simulation.analyses
    results = Results(simulation.model, solution, analyses.subsampler)
    analyse(analyses, results, simulation.output)
end

function results_only(jl_filename::AbstractString)
    include(jl_filename)
    @load "parameters.jld2" simulation
    solution = solve(simulation)
    return (simulation, Results(simulation.model, solution, simulation.analyses.subsampler))
end

export simulate, Simulation, write_params, Solver,
    time_span, solve

end