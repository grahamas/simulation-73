module Modeling

using CalculatedParameters

abstract type Model{T,N,P} <: Parameter{T} end

function space_arr(model::Model)::AbstractArray
    Calculated(model.space).value
end

function initial_value(model::Model{T,N,P}) where {T,N,P}
	space_zeros = zeros(model.space)
end

export Model, initial_value, space_arr

end