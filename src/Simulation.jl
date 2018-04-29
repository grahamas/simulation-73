# * Load Modules
using TensorOperations

using DifferentialEquations

import Base.Dates

# * Simulation function definitions
# ** Core differential functions
#=
Two versions of the differential equation function, both implementing the same
equation described in the Equations section. One is general, putatively
operating on any number of populations (though must still be 1D!). It is much
less efficient in the 2-pop case than the latter, which takes a matrix instead
of a tensor as the interaction term.
=#
function WilsonCowan73!(dA::SpaceState1D{ValueT},A::SpaceState1D{ValueT},
               p::WilsonCowan73Params{InteractionTensor{ValueT},ExpandedParam{ValueT}},
               t::TimeT) where {ValueT <: Real, TimeT <: Real}
    # Use dA as intermediate variable for tensor op since it is preallocated
    @tensor dA[x_tgt, pop_tgt] = p.W[x_tgt, pop_tgt, x_src, pop_src] * A[x_src, pop_src]
    dA .= (-p.α .* A + p.β .* (1 .- A) .* p.nonlinearity_fn(dA + p.stimulus_fn(t))) ./ p.τ
end

doc"""
This simulates a single timestep of the 1D spatial equation
```math
\begin{align}
\tau_E \partial_t E(x,t) &= -\alpha_E E(x,t) + \beta_E (1 - E(x,t)) \cS_E \left( W_{EE}(X) \conv E(x,t) + W_{EI}(X) \conv I(x,t) + P_E(x,t)\right)\\
\tau_I \partial_t I(x,t) &= -\alpha_I I(x,t) + \beta_I (1 - I(x,t)) \cS_I \left( W_{IE}(X) \conv E(x,t) +  W_{II}(X) \conv I(x,t) + P_I(x,t)\right)
\end{align}
```
"""
function WilsonCowan73!(dA::SpaceState1DFlat{ValueT},A::SpaceState1DFlat{ValueT},
                        p::WilsonCowan73Params{Interaction1DFlat{ValueT},ExpandedParamFlat{ValueT}},
                        t::TimeT) where {ValueT <: Real, TimeT<: Real}
    # Use dA as intermediate variable for tensor op since it is preallocated
    # println(size(p.β .* (1 .- A) .* p.nonlinearity_fn(p.W*A + p.stimulus_fn(t))))
    dA .= (-p.α .* A + p.β .* (1 .- A) .* p.nonlinearity_fn(p.W*A + p.stimulus_fn(t))) ./ p.τ
end

# * Solver function definition

function solve_WilsonCowan73(; model=nothing, solver=nothing, other...)
    # The solver parameters are modified, but are used elsewhere so need to be copied
    solver_params = deepcopy(Dict{Any,Any}(solver))
    model_params = WilsonCowan73Params(model)

    u0 = zeros(model_params.mesh)
    tspan = (0.0, pop!(solver_params, :T))
    prob::DEProblem = ODEProblem(WilsonCowan73!, u0, tspan, model_params)

    if :dt in keys(solver_params)
        solver_params[:alg] = Euler()
        solver_params[:adaptive] = false
    elseif :stiff in keys(solver_params)
        if pop!(solver_params, :stiff) > 0
            solver_params[:alg_hints] = [:stiff]
        end
    end

    if :stiff in keys(solver_params)
        error("Incompatible solver parameters.")
    end

    soln::DESolution = solve(prob; solver_params...)

    return soln
end

export solve_WilsonCowan73
