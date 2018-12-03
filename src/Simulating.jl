module Simulating

using Modeling
using Analysis
using Records
using Parameters
using JLD2
import DifferentialEquations: DESolution, solve


@with_kw struct Solver
    T::Float64
    params::Dict
end
time_span(solver::Solver) = (0.0, solver.T)

@with_kw struct Simulation{M<:Model}
    model::M
    solver::Solver
    analyses::Analyses{M}
    output::Output
end

Modeling.initial_value(sim::Simulation) = initial_value(sim.model)
time_span(sim::Simulation) = time_span(sim.solver)
solver_params(sim::Simulation) = sim.solver.params

function write_params(sim::Simulation)
	write_object(sim.output, "parameters.jld2", "sim", sim)
end

function solve(simulation::Simulation)
    problem = generate_problem(simulation)
    params = solver_params(simulation)
    soln = solve(problem; params...)
    return soln
end

function simulate(jl_filename::AbstractString)
	include(jl_filename)
	@load "parameters.jld2" simulation
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

export simulate, Simulation, write_params, Solver,
    time_span, solve

end