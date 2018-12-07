module Modeling

using CalculatedParameters

abstract type Model{T} <: Parameter{T} end

function get_space(model::Model)::AbstractArray
    Calculated(model.space).value
end

initial_value(model::Model) = repeat(zeros(model.space),
                                     outer=([1 for x in 1:ndims(zeros(model.space))]..., length(model.pop_names)))


export Model, initial_value, get_space

end