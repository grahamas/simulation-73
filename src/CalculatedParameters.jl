module CalculatedParameters

using Logging

# * Parameter
abstract type Variable{T<:Real} end

"A parameter to vary, without bounds."
struct UnboundedVariable{T} <: Variable{T}
    value::T
end

struct BoundedVariable{T} <: Variable{T}
    value::T
    bounds::Tuple{T,T}
end
import Base: zero
zero(::Type{Union{BoundedVariable{T}, T}}) where T <: Number = zero(T)
default_value(var::V) where V <: Variable = var.value
bounds(var::V) where V <: Variable = var.bounds

abstract type Parameter{T<:Union{S,Variable{S}} where S<:Real} end

"Map a type over kwargs"
function pops(t::Type{T}; kwargs...)::Array{T} where T
    syms = keys(kwargs)
    args = zip(values(kwargs)...)
    new_kwargs = (zip(syms, arg) for arg in args)
    return map((single_kwargs) -> t(;single_kwargs...), new_kwargs)
end

abstract type CalculatedParam{S <: Parameter} end
function get_value(cp::CalculatedParam)
    cp.value
end

function update(olds::AA, news::AA) where {AA <: AbstractArray}
    AA[update(pair...) for pair in zip(olds, news)]
end

function Calculated(a::T) where T
	error("Calculated type undefined for type $T")
end

export Variable, UnboundedVariable, BoundedVariable
export zero, default_value, bounds
export Parameter, CalculatedParam, pops, Calculated
export update
export get_value

end