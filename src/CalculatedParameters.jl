module CalculatedParameters

# * Parameter
abstract type Variable{T<:Real} end

struct UnboundedVariable{T} <: Variable{T}
    doc"A parameter to vary, without bounds."
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

function pops(t; kwargs...)
    doc"Map a type over kwargs"
    syms, args = zip(values(kwargs)...)
    single_args = zip(args...)
    new_kwargs = (zip(syms, arg) for arg in single_args)
    return map((single_kwargs) -> t(;single_kwargs...), new_kwargs)
end

abstract type CalculatedParam{S <: Parameter} end

function update!(olds::Array, news::Array)
    any(update!(pair...) for pair in zip(olds, news))
end

function Calculated(a::T) where T
	error("Calculated type undefined for type $T")
end

export Variable, UnboundedVariable, BoundedVariable
export zero, default_value, bounds
export Parameter, CalculatedParam, pops, Calculated
export update!

end