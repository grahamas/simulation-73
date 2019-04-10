

"A model with variable parameters, and a target."
struct ParameterSearch{T, M<: AbstractModel{<:MaybeVariable{T}}, S<:Solver{T}}
    model::M
    solver::S
    target::AbstractTarget
    result::Union{BlackBoxOptim.OptimizationResults, Optim.OptimizationResults}
    result_simulation::Simulation{T,<:AbstractModel{T},S}
end

minimizer(res::BlackBoxOptim.OptimizationResults) = best_candidate(res)
minimizer(res::Optim.OptimizationResults) = Optim.minimizer(res)

function ParameterSearch(;varying_model::M=nothing, solver::S=nothing,
        target::AbstractTarget=nothing, kwargs...) where {T,M<:AbstractModel{<:MaybeVariable{T}},S<:Solver{T}}
    initial_p, variable_dxs, p_bounds = init_variables(varying_model)
    result = run_search(varying_model, variable_dxs, solver, target, initial_p, p_bounds; kwargs...)
    result_model = model_from_p(varying_model, variable_dxs, minimizer(result))
    @show result_model
    result_simulation = Simulation(; model = result_model, solver=solver)
    @info "Left simulation result"
    ParameterSearch{T,M,S}(varying_model, solver, target, result, result_simulation)
end

result_simulation(p_search::ParameterSearch) = p_search.result_simulation

"""Takes model with variable parameters,
and returns default variable values and indices to those variables."""
function init_variables(variable_model::M) where {T,M <: AbstractModel{<:MaybeVariable{T}}}
    deconstructed = var_deconstruct(variable_model)
    initial_p, variable_dxs, p_bounds = init_variables(deconstructed, T)
    return initial_p, variable_dxs, p_bounds
end

"Initialize all variables to their default value."
function init_variables(deconstructed::Tuple{Type,<:AbstractArray}, T::Type)
    variable_dxs = []
    initial_p = T[]
    p_bounds = Tuple{T,T}[]
    for dx in CartesianIndices(size(deconstructed[2]))
        (typ, val) = deconstructed[2][dx]
        if typ <: AbstractVariable
            push!(variable_dxs, [dx])
            push!(initial_p, default_value(val))
            push!(p_bounds, bounds(val))
        elseif val isa AbstractArray
            ps, dxs, bds = init_variables((typ, val), T)
            @assert length(ps) == length(dxs)
            push!(variable_dxs, [vcat(dx, x) for x in dxs]...)
            push!(initial_p, ps...)
            push!(p_bounds, bds...)
        end
    end
    return initial_p, variable_dxs, p_bounds
end

function time_span(p_search::ParameterSearch)
    time_span(p_search.solver)
end



"""
    model_from_p(p_seach::ParameterSearch, p)

Deconstructs the "variable model" stored by p_search into a flat representation.
Indexes into the flat representation using p_search's variable_map, which was
created to relate the locations of varying parameters in p_search's variable
model to the parameters in the parameter vector p.
"""
function model_from_p(model, variable_map, new_p) # Does the model need to be varying? No?
    model_deconstruction = deconstruct(model)
    for (dx, val) in enumerate(new_p)
        target_dxs = variable_map[dx]
        set_deep_dx!(model_deconstruction, target_dxs, (typeof(val), val))
    end
    return reconstruct(model_deconstruction)
end

# * Run the search

make_problem_generator() = error("undefined.")
export make_problem_generator

function write_params(p_search::ParameterSearch)
    write_object(p_search.output, "search_parameters.jld2", "p_search", p_search)
end

function write_results(p_search::ParameterSearch, results)
    write_object(p_search.output, "results.jld2", "results", results)
end

function _build_loss_objective(initial_problem, solver::Solver{T,Euler}, space, loss_fn, prob_generator) where T
    build_loss_objective(initial_problem, Euler(), loss_fn;
        prob_generator=prob_generator, dt=solver.simulated_dt,
        save_at=save_dt(solver), save_idxs=save_idxs(solver, space))
end

function _build_loss_objective(initial_problem, solver::Solver{T,Nothing}, space, loss_fn, prob_generator) where T
    build_loss_objective(initial_problem, Tsit5(), loss_fn;
        prob_generator=prob_generator,
        save_at=save_dt(solver),
        timeseries_steps=solver.time_save_every,
        save_idxs=save_idxs(solver, space),
        alg_hints=[solver.stiffness])
end

function run_search(varying_model, variable_map, solver, target, initial_p, p_bounds::Array{<:Tuple}; MaxSteps=11e3)
    @show MaxSteps
    initial_model = model_from_p(varying_model, variable_map, initial_p)
    problem_generator = make_problem_generator(initial_model, solver, variable_map)
    initial_problem = problem_generator(nothing, initial_p)
    loss_fn = target_loss(target, initial_model, solver)
    loss_obj = _build_loss_objective(initial_problem, solver, initial_model.space, loss_fn,
                                                problem_generator)
    result = bboptimize(loss_obj; NumDimensions=length(initial_p),
        MaxSteps=MaxSteps, SearchRange=p_bounds)
    return result
end

function run_search_optim(varying_model, variable_map, solver, target, initial_p, p_bounds::Array{<:Tuple}; MaxSteps=11e3)
    @show MaxSteps
    initial_model = model_from_p(varying_model, variable_map, initial_p)
    problem_generator = make_problem_generator(initial_model, solver, variable_map)
    initial_problem = problem_generator(nothing, initial_p)
    loss_fn = target_loss(target, initial_model, solver)
    loss_obj = _build_loss_objective(initial_problem, solver, initial_model.space, loss_fn,
                                                problem_generator)
    lower = [bounds[1] for bounds in p_bounds]
    upper = [bounds[2] for bounds in p_bounds]
    result = optimize(loss_obj, lower, upper, initial_p, Fminbox(ParticleSwarm()),
                        Optim.Options(
                         iterations = 10,
                         show_every = 1))
    return result
end

function run_search(jl_filename::AbstractString)
    include(jl_filename)
    filecopy(p_search.output, jl_filename, basename(jl_filename))
    return p_search
end

update_from_p!(args...) = error("Please import and implement.")
make_calculated_function(args...) = error("Please import and implement.")

function make_problem_generator(model::M, solver::SV, variable_map) where {T,M<:AbstractModel{T},SV<:Solver{T}}
    tspan = time_span(solver)
    u0 = initial_value(model)

    calculated_model = Calculated(model)

    function problem_generator(prob, new_p::Array{T})
        update_from_p!(calculated_model, new_p, model, variable_map)
        calculated_function = make_calculated_function(calculated_model)
        ode_fn = convert(ODEFunction{true}, calculated_function)
        return ODEProblem(ode_fn, u0, tspan, new_p)
    end

    return problem_generator
end
