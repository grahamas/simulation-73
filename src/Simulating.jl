module Simulating

using Modeling
using Analysis
using Records
using Parameters


@with_kw struct Solver
    T::Float64
    params::Dict
end
time_span(solver::Solver) = (0.0, solver.T)

function solve(problem, solver::Solver)
	solve(problem; solver.params...)
end


@with_kw struct Simulation{M<:Model}
    model::M
    solver::Solver
    analyses::Analyses{M}
    output::Output
end

Modeling.initial_value(sim::Simulation) = initial_value(sim.model)
time_span(sim::Simulation) = time_span(sim.solver)

function write_params(sim::Simulation)
	write_object(sim.output, "parameters.jld2", "sim", sim)
end


export Simulation, write_params, Solver, time_span, solve

end