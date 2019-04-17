
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

"""
    deconstruct(ad::AbstractDeconstructor, param)

Deconstruct a parameter into nested tuples of the form `(type,value)`.

```
    Parameter{Varying{T}} |> deconstruct |> reconstruct # isa Parameter{Varying{T}}
```
"""
function deconstruct(m::M, ad::AbstractDeconstructor=Deconstructor()) where {M <: AbstractParameter}
    deconstruction = Tuple{Type,Any}[]
    for i_field in 1:nfields(m)
        substruct = ad(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return (typeof(m), deconstruction)
end

"""
	deconstruct(..., DeconstructorFixingVariables())

Deconstructs as normal (i.e. [Deconstructor](@ref)) except on [AbstractVariable](@ref) parameters, which are fixed to their default value.
"""
struct DeconstructorFixingVariables <: AbstractDeconstructor end

function deconstruct(val::AbstractVariable, ad::AbstractDeconstructor=Deconstructor())
    return (typeof(val), val)
end
function deconstruct(v::AbstractVariable, deconstructor::DeconstructorFixingVariables)
    val = default_value(v)
    return (typeof(val), val)
end


function base_type(::Type{P}) where {P <: AbstractParameter}
    BTs = map(base_type, P.parameters)
    return (P.name.wrapper){BTs...}
end

function base_type(::Type{V}) where {T, V<: Union{AbstractVariable{T}, T}}
    return T
end
