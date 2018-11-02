module Modeling

# * Load Modules
using DifferentialEquations
import DifferentialEquations: solve, DESolution#, DEProblem
using Records

import Dates

using Parameters
using CalculatedParameters
# * Simulation object

abstract type Model{T} <: Parameter{T} end
function space(model::Model)::AbstractArray
    Calculated(model.space).value
end

initial_value(model::Model) = repeat(zeros(model.space),
                                     outer=([1 for x in 1:ndims(zeros(model.space))]..., length(model.pop_names)))


export Model
export initial_value

end