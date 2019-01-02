module WCMNonlinearity
# * Types

using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update!

abstract type Nonlinearity{T} <: Parameter{T} end

function update!(calc_nln::Array{CL,1}, new_nln::Array{L,1}) where {T, L <: Nonlinearity{T}, CL<:CalculatedParam{L}}
    for i in 1:length(calc_nln)
        if calc_nln[i].nonlinearity != new_nln[i]
            calc_nln[i] = Calculated(new_nln[i])
        end
    end
end

@with_kw struct SigmoidNonlinearity{T} <: Nonlinearity{T}
    a::T
    θ::T
end

mutable struct CalculatedSigmoidNonlinearity{T} <: CalculatedParam{SigmoidNonlinearity{T}}
    nonlinearity::SigmoidNonlinearity{T}
    a::T
    θ::T # This is nonsense, but fits with Calculated pattern
end

function Calculated(sigmoid::SigmoidNonlinearity{T}) where T
    CalculatedSigmoidNonlinearity{T}(sigmoid, sigmoid.a, sigmoid.θ)
end

function nonlinearity!(output::AT, csn::CalculatedSigmoidNonlinearity{T}) where {T, AT<:AbstractArray{T}}
    output .= rectified_sigmoid_fn.(output,csn.a,csn.θ)
end

CalculatedParameters.get_value(csn::CalculatedSigmoidNonlinearity{T}) where T = csn

function sech2_fn(x, a, θ)
    return @. 1 - tanh(a * (x - θ))^2
end

"""
The sigmoid function is defined

```math
\\begin{align}
\\mathcal{S}(x) = \\frac{1}{1 + \\exp(-a(x - θ))}
\\end{align}
```
where ``a`` describes the slope's steepness and ``θ`` describes translation of the slope's center away from zero.

This is "simple" because in practice we use the rectified sigmoid.
"""
function simple_sigmoid_fn(x, a, theta)
    return @. (1.0 / (1 + exp(-a * (x - theta))))
end

"""
A rectified version of `simple_sigmoid_fn`.

In practice, we use rectified sigmoid functions because firing rates cannot be negative.

TODO: Rename to rectified_sigmoid_fn.
"""
function rectified_sigmoid_fn(x, a, theta)
    return max.(0, simple_sigmoid_fn(x, a, theta) .- simple_sigmoid_fn(0, a, theta))
end

function sigmoid_diff_fn(x, a, θ, width)
    return max.(0,sigmoid_fn(x, a, θ) - sigmoid_fn(x, a, θ + width))
end

function neg_domain_sigmoid_diff_fn(input, a, θ, width)
    return max.(0,simple_sigmoid_fn(input, a, θ) - simple_sigmoid_fn(input, a, θ + width))
end

export Nonlinearity, SigmoidNonlinearity, CalculatedSigmoidNonlinearity, Calculated, update!, nonlinearity!

# * end
end
