module Exploration

using Parameters
using Modeling
using Analysis
using Records
using CalculatedParameters
using Targets
using DifferentialEquations
using BlackBoxOptim
import Records: required_modules

import BlackBoxOptim: OptimizationResults
required_modules(::Type{OptimizationResults}) = [BlackBoxOptim]

doc"A model with variable parameters, and a target."
struct ParameterSearch{M<: Model}
    model::M
    solver::Solver
    analyses::Analyses
    output::Output
    target::Target
    initial_p
    variable_map
    p_bounds
end

function required_modules(::Type{ParameterSearch{M}}) where {M <: Model}
    union([Exploration, Records, Analysis, Targets], required_modules(M))
end

function ParameterSearch(;variable_model::M=nothing, solver::Solver=nothing,
                     analyses::Analyses=nothing, output::Output=nothing,
                     target::Target=nothing) where {M<:Model}
    initial_p, variable_dxs, p_bounds = init_variables(variable_model)
    ParameterSearch{M}(variable_model, solver, analyses, output, target,
        initial_p, variable_dxs, p_bounds)
end

export ParameterSearch, required_modules

function init_variables(variable_model::M) where {M <: Model}
    doc"Takes model with variable parameters, and returns default variable values and indices to those variables."
    deconstructed = deconstruct(variable_model)
    initial_p, variable_dxs, p_bounds = init_variables(deconstructed)
    return initial_p, variable_dxs, p_bounds
end

function init_variables(deconstructed::Tuple{Type,<:AbstractArray})
    variable_dxs = []
    initial_p = []
    p_bounds = Tuple{Float64,Float64}[]
    for dx in CartesianRange(size(deconstructed[2]))
        (typ, val) = deconstructed[2][dx]
        if typ <: Variable
            push!(variable_dxs, [dx])
            push!(initial_p, default_value(val))
            push!(p_bounds, bounds(val))
        elseif val isa AbstractArray
            ps, dxs, bds = init_variables((typ, val))
            @assert length(ps) == length(dxs)
            push!(variable_dxs, [vcat(dx, x) for x in dxs]...)
            push!(initial_p, ps...)
            push!(p_bounds, bds...)
        end
    end
    return initial_p, variable_dxs, p_bounds
end

function initial_model(p_search::ParameterSearch)
    model_from_p(p_search, p_search.initial_p)
end

import Modeling: time_span
function time_span(p_search::ParameterSearch)
    time_span(p_search.solver)
end

export initial_model, time_span

function set_deep_dx!(model_deconstruction::Tuple{Type,Array}, dxs, val)
    tmp = model_deconstruction
    for dx in dxs[1:end-1]
        tmp = tmp[2][dx]
    end
    tmp[2][dxs[end]] = val
end

function deconstruct(val::O) where {O <: Union{Real, Variable}}
    return (typeof(val), val)

end

function deconstruct(m::M) where {M <: Parameter}
    deconstruction = Tuple{Type,Any}[]
    for i_field in 1:nfields(m)
        substruct = deconstruct(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return (typeof(m), deconstruction)
end

function deconstruct(arr::Array)
    deconstruction = map(deconstruct, arr)
    return (base_type(typeof(arr)), deconstruction)
end

function base_type(::Type{P}) where {P <: Parameter}
    BTs = map(base_type, P.parameters)
    return (P.name.wrapper){BTs...}
end

function base_type(::Type{T}) where T <: Real
    return T
end

function base_type(::Type{V}) where {T, V<: Union{Variable{T}, T}}
    return T
end

function base_type(::Type{Array{T,N}}) where {T, N}
    BT = base_type(T)
    return Array{BT,N}
end

function reconstruct(tup::Tuple{Type,Real})
    return tup[2]
end

function reconstruct(tup::Tuple{Type,Array})
    typ, arr = tup
    base_typ = base_type(typ)
    if typ <: Array
        return base_typ(reconstruct.(arr))
    else
        return base_typ(reconstruct.(arr)...)
    end
end

function model_from_p(p_search::ParameterSearch, p)
    var_model = p_search.model
    variable_map = p_search.variable_map
    model_deconstruction = deconstruct(var_model)
    for (dx, val) in enumerate(p)
        target_dxs = variable_map[dx]
        set_deep_dx!(model_deconstruction, target_dxs, (typeof(val), val))
    end
    return reconstruct(model_deconstruction)
end

function simulation_from_p(p_search::ParameterSearch, p)
    model = model_from_p(p_search, p)
    return Simulation(model, p_search.solver,
        p_search.analyses, p_search.output)
end
# * Run the search

function loss(factory::TargetFactory, model::Model)
    target_data = factory(model)
    weights = zeros(target_data)
    weights[Calculated(model.space).value .>= 0,factory.target_pop,:] .= 1
    L2Loss(factory.timepoints, target_data; data_weight=weights)
end

function loss(fn::LossFunction, model::Model)
    calc_space = Calculated(model.space).value
    return (soln) -> fn(soln, calc_space)
end

make_problem_generator() = error("undefined.")
export make_problem_generator

function write_params(p_search::ParameterSearch)
    write_object(p_search.output, "search_parameters.jld2", "p_search", p_search)
end

function write_results(p_search::ParameterSearch, results)
    write_object(p_search.output, "results.jld2", "results", results)
end

function run_search(p_search::ParameterSearch)
    write_params(p_search)
    initial_problem, problem_generator = make_problem_generator(p_search)
    loss_fn = loss(p_search.target, initial_model(p_search))
    loss_obj = build_loss_objective(initial_problem, Tsit5(), loss_fn;
                                    prob_generator=problem_generator)
    result = bboptimize(loss_obj; NumDimensions=length(p_search.initial_p),
        MaxSteps=1, SearchRange=p_search.p_bounds)
    write_results(p_search, result)
    result_sim = simulation_from_p(p_search, best_candidate(result))
    result_problem = problem_generator(nothing, best_candidate(result))
    result_sol = solve(result_problem, p_search.solver)
    analyse(result_sim, result_sol)
    return result
end

import DifferentialEquations: solve

function solve(simulation::Simulation)
    initial_problem, problem_generator = make_problem_generator(simulation)
    params = solver_params(simulation)
    soln = solve(initial_problem; params...)
    return soln
end

export solve, run_search, model_from_p

end
