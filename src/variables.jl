


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


# const UV = UnboundedVariable
# const BV = BoundedVariable
# const Varying{T} = Union{T,BV{T}}

# export Varying



"Map a type over kwargs"
function pops(t::Type{T}; kwargs...)::Array{T} where T
    syms = keys(kwargs)
    args = zip(values(kwargs)...)
    new_kwargs = (zip(syms, arg) for arg in args)
    return map((single_kwargs) -> t(;single_kwargs...), new_kwargs)
end
