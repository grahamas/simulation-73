module CalculatedParameter

abstract type Parameter{T} where {T <: Real} end

function pops(t; kwargs...)
    doc"Map a type over kwargs"
    return map((v) -> t(Dict(zip(keys(kwargs), v))...), zip(values(kwargs)...))
end

abstract type Calculated{T} where {T <: Parameter} end

function update!(olds::Array{<:Calculated{T}}, news::Array{T}) where {T <: Parameter}
    any(update!(pair...) for pair in zip(olds, news))
end

export Parameter, Calculated, pops

end
