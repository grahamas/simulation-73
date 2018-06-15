
module ParameterSearch
# * Parameter
abstract type Variable end

struct UnboundedVariable{T<:Real}
    value::T
end

function variable_p(x)
    dxs = Int[]
    ps = []
    for dx in 1:nfields(x)
        val = variable_p(getfield(x,dx))
        if val != nothing
            push!(dxs, dx)
            push!(ps, val)
        end
    end
    return dxs, ps
end

function variable_p(x::Variable)
    return x.value
end

function variable_p(x::Number)
    return nothing
end

# * Parameter search type
struct ParameterSearch{M <: Model}
    model::M
    solver::Solver
    analyses::Analyses
    output::Output
    target::Target
    variable_map
end

function fix_variables(variable_model::M) where {T<:Real,M <: Model{Union{Variable,T}}}
    pass
end

function ParameterSearch(variable_model::Model{Union{Variable, T<:Real}}, solver::Solver,
                         analyses::Analyses, output::Output)
    variable_map, model = fix_variables(variable_model)
    ParameterSearch{M{T}}(model, solver, analyses, output, variable_map)
end

# ** Helpers to modify parameters

function set_deep_dx!(array, dxs, val)
    tmp = array
    for dx in dxs[1:end-1]
        tmp = tmp[dx]
    end
    tmp[dxs[end]] = val
end

function deconstruct(val::Real)
    val
end

function deconstruct(m)
    deconstruction = []
    for i_field in 1:nfields(m)
        substruct = deconstruct(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return deconstruction
end

function reconstruct(m, val::Real)
    val
end

function reconstruct(m, arr)
    typeof(m)([reconstruct(getfield(m, i), arr[i]) for i in 1:nfields(m)]...)
end

function model_from_p(p_search::ParameterSearch, p)
    model = p_search.model
    model_array = deconstruct(model)
    variable_map = p_search.variable_map
    for (dx, val) in enumerate(p)
        target_dxs = variable_map[dx]
        set_deep_dx!(model_array, target_dxs, val)
    end
    return reconstruct(model, model_array)
end

# * Run the search

function run(search::ParameterSearch)
    initial_problem, problem_generator = make_problem_generator(search)
    loss_fn = loss(target, search.model)
    loss_obj = build_loss_objective(initial_problem, Tsit5(), loss_fn;
                                    prob_generator=problem_generator)
end

function solve(simulation::WC73Simulation)
    initial_problem, problem_generator = make_problem_generator(simulation)
    params = solver_params(simulation)
    soln = solve(initial_problem; params...)
    return soln
end

export solve

end
