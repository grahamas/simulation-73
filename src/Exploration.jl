
module Exploration

using Parameters
using Simulation

# * Parameter
abstract type Variable{T<:Real} end

struct UnboundedVariable{T<:Real} <: Variable{T}
    doc"A parameter to vary, without bounds."
    value::T
end

struct BoundedVariable{T<:Real} <: Variable{T}
    value::T
    bounds::Tuple{T,T}
end

doc"""Takes model to make true target."""
@with_kw struct TargetFactory{T<:Real}
    factory::Function
    timepoints::Array{T,1}
end


doc"A model with variable parameters, and a target."
struct ParameterSearch{M<: Model{Union{Variable{T<:Real},T}}}
    model::M
    solver::Solver
    analyses::Analyses
    output::Output
    target::TargetFactory
    initial_p
    variable_map
end

function ParameterSearch(;variable_model::M, solver::Solver,
                     analyses::Analyses, output::Output,
                     target_factory::TargetFactory) where {M<:Model}
    initial_p, variable_dxs = init_variables(variable_model)
    model = model_from_p(variable_model, initial_p)
    ParameterSearch{M}(variable_model, solver, analyses, output, target, initial_p, variable_dxs)
end

function init_variables(variable_model::M) where {T<:Real,M <: Model{Union{Variable{T},T}}}
    doc"Takes model with variable parameters, and returns default variable values and indices to those variables."
    deconstructed = deconstruct(variable_model)
    initial_p, variable_dxs = init_variables(deconstructed)
    return intial_p, variable_dxs
end

function init_variables(deconstructed::Array{Tuple{Type,Array}})
    variable_dxs = []
    initial_p = []
    for (dx, (typ, val)) in deconstructed
        if typ <: Variable
            p = default_value(typ(val))
            push(variable_dxs, [dx])
            push(initial_p, p)
        else
            ps, dxs = (typ, init_variables(val))
            @assert length(ps) == length(dxs)
            push!(variable_dxs, [vcat(dx, x) for x in dxs]...)
            push!(initial_p, ps...)
        end
    end
    return initial_p, variable_dxs
end

function init_variables(arr::Array)
    return arr, []
end

function initial_model(p_search::ParameterSearch)
    model_from_p(p_search, p_search.initial_p)
end

function time_span(p_search::ParameterSearch)
    time_span(p_search.solver)
end

function set_deep_dx!(model_deconstruction::Tuple{Type,Array}, dxs, val)
    tmp = model_deconstruction
    for dx in dxs[1:end-1]
        tmp = tmp[2][dx]
    end
    tmp[dxs[end]] = val
end

function deconstruct(val::Real)
    return val
end

function deconstruct(m)
    deconstruction = Tuple{Type,Array}[]
    for i_field in 1:nfields(m)
        substruct = deconstruct(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return (typeof(m), deconstruction)
end

function deconstruct(arr::Array)
    deconstruction = []
    for (dx, val) in enumerate(arr)
        substruct = deconstruct(val)
        push!(deconstruction, substruct)
    end
    return (typeof(arr), deconstruction)
end

function reconstruct(val::Real)
    return val
end

function reconstruct(tup::Tuple{Type,Array})
    typ, arr = tup
    if typ <: Array
        return typ(reconstruct.(arr))
    else
        return typ(reconstruct.(arr)...)
    end
end

function model_from_p(var_model::Model,p)
    model_deconstruction = deconstruct(var_model)
    variable_map = p_search.variable_map
    for (dx, val) in enumerate(p)
        target_dxs = variable_map[dx]
        set_deep_dx!(model_deconstruction, target_dxs, val)
    end
    return reconstruct(model_deconstruction)
end

function model_from_p(p_search::ParameterSearch, p)
    model_from_p(p_search.model, p)
end

function simulation_from_p(p_search::ParameterSearch, p)
    model = model_from_p(p_search, p)
    return Simulation(model, p_search.solver,
        p_search.analyses, p_search.output)
end
# * Run the search

function loss(factory::TargetFactory, model::Model)
    timepoints = factory.timepoints
    target_fn = factory.factory(model)
    target_data = target_fn.(timepoints)
    L2Loss(timepoints, target_data)
end

function run(p_search::ParameterSearch)
    initial_problem, problem_generator = make_problem_generator(p_search)
    loss_fn = loss(p_search.target, p_search.model)
    loss_obj = build_loss_objective(initial_problem, Tsit5(), loss_fn;
                                    prob_generator=problem_generator)
    result = bboptimize(loss_obj)
    result_problem = problem_generator(nothing, best_candidate(result))
    result_sol = solve(result_problem, p_search.solver)
    analyse_WilsonCowan73_solution(result_problem, p_search.analyses)
    return result
end

function solve(simulation::WC73Simulation)
    initial_problem, problem_generator = make_problem_generator(simulation)
    params = solver_params(simulation)
    soln = solve(initial_problem; params...)
    return soln
end

export solve

end
