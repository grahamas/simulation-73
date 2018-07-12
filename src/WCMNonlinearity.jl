module WCMNonlinearity
# * Types

using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update!

abstract type Nonlinearity{T} <: Parameter{T} end

@with_kw struct SigmoidNonlinearity{T} <: Nonlinearity{T}
    a::T
    θ::T
end

function calculate(sn::SigmoidNonlinearity)
    make_sigmoid_fn(sn.a, sn.θ)
end

mutable struct CalculatedSigmoidNonlinearity{T} <: CalculatedParam{SigmoidNonlinearity{T}}
    sigmoid::SigmoidNonlinearity{T}
    value::Function
    CalculatedSigmoidNonlinearity{T}(s::SigmoidNonlinearity{T}) where T<:Number = new(s, calculate(s))
end

function Calculated(sigmoid::SigmoidNonlinearity{T}) where T
    CalculatedSigmoidNonlinearity{T}(sigmoid)
end

function update!(csn::CalculatedSigmoidNonlinearity, sn::SigmoidNonlinearity)
    if csn.sigmoid == sn
        return false
    else
        csn.sigmoid = sn
        csn.value = calculate(sn)
    end
end

# * Sech2
function make_sech2_fn(; a=error("Missing arg"), θ=error("Missing arg"))
    return (x) -> max.(0,sech2_fn(x, a, θ))
end

function sech2_fn(x, a, θ)
    return @. 1 - tanh(a * (x - θ))^2
end

# * Sigmoid functions
function make_sigmoid_fn(a, θ)
    return (x) -> sigmoid_fn(x, a, θ)
end
doc"""
The sigmoid function is defined
```math
\begin{align}
\mathcal{S}(x) = \frac{1}{1 + \exp(-a(x - \theta))}
\end{align}
```
where $a$ describes the slope's steepness and $\theta$ describes translation of the slope's center away from zero.

This is "simple" because in practice we use the rectified sigmoid.
"""
function simple_sigmoid_fn(x, a, theta)
    return @. (1.0 / (1 + exp(-a * (x - theta))))
end

doc"""
A rectified version of `simple_sigmoid_fn`.

In practice, we use rectified sigmoid functions because firing rates cannot be negative.

TODO: Rename to rectified_sigmoid_fn.
"""
function sigmoid_fn(x, a, theta)
    return max.(0, simple_sigmoid_fn(x, a, theta) .- simple_sigmoid_fn(0, a, theta))
end

# * Difference of sigmoids functions
apply(fn, x) = fn(x)
function make_sigmoid_diff_fn(; a=nothing, θ=nothing, width=nothing)
    if size(a) == ()
        return make_sigmoid_diff_fn(a, θ, width)
    end
    arg_list = collect(zip(a, θ, width))
    fn_list = arg_list .|> (args) -> make_sigmoid_diff_fn(args...)
    return (xs) -> apply.(fn_list, xs)
end

function make_sigmoid_diff_fn(a::T, θ::T, width::T) where {T <: Real}
    unscaled(y) = sigmoid_diff_fn(y, a, θ, width)  # Peak is not always 1
    range = (θ-(1.0 ./ a)):0.001:(θ+(1.0 ./ a)+width)
    maxes = maximum(unscaled.(range), 1)
    return (x) -> unscaled(x) ./ maxes[1]
end

function make_neg_domain_sigmoid_diff_fn(a::T, θ::T, width::T) where {T <: Real}
    unscaled(y) = neg_domain_sigmoid_diff_fn(y, a, θ, width)  # Peak is not always 1
    range = (θ-(1.0 ./ a)):0.001:(θ+(1.0 ./ a)+width)
    maxes = maximum(unscaled.(range), 1)
    return (x) -> unscaled(x) / maxes[1]
end

function sigmoid_diff_fn(x, a, θ, width)
    return max.(0,sigmoid_fn(x, a, θ) - sigmoid_fn(x, a, θ + width))
end

function neg_domain_sigmoid_diff_fn(input, a, θ, width)
    return max.(0,simple_sigmoid_fn(input, a, θ) - simple_sigmoid_fn(input, a, θ + width))
end

export Nonlinearity, SigmoidNonlinearity, CalculatedSigmoidNonlinearity, Calculated, update!

# * end
end
