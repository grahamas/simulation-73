
abstract type AbstractParameter{T} end
DrWatson.default_allowed(c::AbstractParameter) = (Real, String, Symbol, AbstractParameter)


