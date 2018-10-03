module Modeling

# * Load Modules
using DifferentialEquations
import DifferentialEquations: solve, DESolution#, DEProblem
using Analysis
import Analysis: analyse
using Records

import Dates

using Parameters
using CalculatedParameters
# * Simulation object

abstract type Model{T} <: Parameter{T} end

@with_kw struct Solver
    T::Float64
    params::Dict
end

@with_kw struct Simulation{M<:Model}
    model::M
    solver::Solver
    analyses::Analyses
    output::Output
end

time_span(solver::Solver) = (0.0, solver.T)
time_span(sim::Simulation) = time_span(sim.solver)

initial_value(model::Model) = repeat(zeros(model.space),
                                     outer=([1 for x in 1:ndims(zeros(model.space))]..., length(model.α)))
initial_value(sim::Simulation) = initial_value(sim.model)

function solve(problem, solver::Solver)
	solve(problem; solver.params...)
end

function write_params(sim::Simulation)
	write_object(sim.output, "parameters.jld2", "sim", sim)
end

function analyse(sim::Simulation, solution::DESolution)
	write_params(sim)
	analyse(solution, sim.output, sim.model; sim.analyses...)
end

export Model, Solver, Simulation
export time_span, initial_value, solve, analyse, write_params

end