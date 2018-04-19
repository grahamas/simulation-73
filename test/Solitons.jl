#=
A module for generating solitons as ground truth for testing.
=#

module Solitons

simple_soliton(x,t,c,x_0) = c/2 * sech(sqrt(c)/2 * (x - c * t - x_0)) ^ 2
function make_simple_soliton(c, x_0)
    return (x,t) -> c/2 * sech(sqrt(c)/2 * (x - c * t - x_0)) ^ 2
end

generic_soliton(x,t,A,b,c,x_0) = A * sech(b * x - c * t - x_0) ^ 2
function make_generic_soliton(A, b, c, x_0)
    return (x,t) ->  A * sech(b * x - c * t - x_0) ^ 2
end

exponential_decay(t,A,σ) = A * exp(t / σ)
function make_exponential_decay(A, σ)
    return (t) -> A * exp(t / σ)
end

function time_dependent_soliton(x,t,A,C,s_A,s_C,b,x_0)
    return exponential_decay(t,A,s_A) * sech(b*x - exponential_decay(t,C,s_C)*t - x_0)^2
end
function make_time_dependent_soliton(A,C,s_A,s_C,b,x_0)
    edA = make_exponential_decay(A,s_A)
    edC = make_exponential_decay(C,s_C)
    return (x,t) -> edA(t) * sech(b*x - edC(t)*t - x_0)^2
end

export simple_soliton, make_simple_soliton
export generic_soliton, make_generic_soliton
export exponential_decay, make_exponential_decay
export time_dependent_soliton, make_time_dependent_soliton

end
