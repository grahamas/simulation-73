
abstract type AbstractVariable{T<:Real} end

"A parameter to vary, without bounds."
struct UnboundedVariable{T} <: AbstractVariable{T}
    value::T
end

struct BoundedVariable{T} <: AbstractVariable{T}
    value::T
    bounds::Tuple{T,T}
end

function BoundedVariable(value::T) where {T <: Real}
	ten_pct = value / 10
	BoundedVariable{T}(value, (value-ten_pct,value+ten_pct))
end

const MaybeVariable{T} = Union{T,AbstractVariable{T}}

Base.zero(::Type{<:MaybeVariable{T}}) where T <: Number = zero(T)
default_value(var::V) where V <: AbstractVariable = var.value
bounds(var::V) where V <: BoundedVariable = var.bounds

abstract type AbstractParameter{T<:MaybeVariable} end

"Map a type over kwargs"
function pops(t::Type{T}; kwargs...)::Array{T} where T
    syms = keys(kwargs)
    args = zip(values(kwargs)...)
    new_kwargs = (zip(syms, arg) for arg in args)
    return map((single_kwargs) -> t(;single_kwargs...), new_kwargs)
end