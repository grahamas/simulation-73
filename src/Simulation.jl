module Simulation

# * Load Modules
using TensorOperations

using DifferentialEquations

using ..Analysis

import Base.Dates
# * Simulation object

abstract type Model end

@with_kw struct Simulation{M<:Model}
    model::M
    solver::Solver
    analyses::Analyses
end

@with_kw struct Solver
    T::Float64
    params::Dict
end

time_span(solver::Solver) = (0.0, solver.T)
time_span(sim::Simulation) = time_span(sim.solver)

initial_value(model::Model) = repeat(zeros(model.space),
                                     outer=([1 for x in 1:ndims(zeros(model.space))]..., model.n_pops))
initial_value(sim::Simulation) = initial_value(sim.model)

end
