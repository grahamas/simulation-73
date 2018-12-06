module Simulating

using Modeling
using Analysis
using Records
using Parameters
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve


@with_kw struct Solver
    T::Float64
    params::Dict
    solution_method::Union{OrdinaryDiffEqAlgorithm,Nothing}=nothing
end
time_span(solver::Solver) = (0.0, solver.T)

@with_kw struct Simulation{M<:Model}
    model::M
    solver::Solver
    analyses::Analyses
    output::AbstractOutput
end

Modeling.initial_value(sim::Simulation) = initial_value(sim.model)
time_span(sim::Simulation) = time_span(sim.solver)
solver_params(sim::Simulation) = sim.solver.params
solution_method(sim::Simulation) = sim.solver.solution_method

function write_params(sim::Simulation)
	write_object(sim.output, "parameters.jld2", "sim", sim)
end

function solve(simulation::Simulation)
    problem = generate_problem(simulation)
    params = solver_params(simulation)
    if solution_method(simulation) != nothing
        soln = solve(problem, solution_method(simulation); params...)
    else
        soln = solve(problem; params...)
    end
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
    return results
    analyse(analyses, results, simulation.output)
end

export simulate, Simulation, write_params, Solver,
    time_span, solve

end