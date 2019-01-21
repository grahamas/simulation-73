module Exploration

using Parameters
using Modeling
using Simulating
using Analysis
using Records
using CalculatedParameters
using Targets
using DifferentialEquations, DiffEqParamEstim
using BlackBoxOptim
import Records: required_modules
using StaticArrays

const UV = UnboundedVariable
const BV = BoundedVariable
const Varying{T} = Union{T,BV{T}}

export Varying

import BlackBoxOptim: OptimizationResults
required_modules(::Type{OptimizationResults}) = [BlackBoxOptim]

"A model with variable parameters, and a target."
struct ParameterSearch{T, M<: Model{Varying{T}}, S<:Solver{T}}
    model::M
    solver::S
    analyses::Analyses
    output::AbstractOutput
    target::Target
    result::OptimizationResults
    simulation_result::Simulation{T,<:Model{T},S}
end

function ParameterSearch(;varying_model::M=nothing, solver::S=nothing,
                     analyses::Analyses=nothing, output::AbstractOutput=nothing,
                     target::Target=nothing, kwargs...) where {T,M<:Model{Varying{T}},S<:Solver{T}}
    initial_p, variable_dxs, p_bounds = init_variables(varying_model)
    result = run_search(varying_model, variable_dxs, solver, target, initial_p, p_bounds; kwargs...)
    model_result = model_from_p(varying_model, variable_dxs, best_candidate(result))
    @show model_result
    simulation_result = Simulation(; model = model_result, solver=solver, analyses=analyses, output=output)
    ParameterSearch{T,M,S}(varying_model, solver, analyses, output, target,
        result, simulation_result)
end

export ParameterSearch

"""Takes model with variable parameters,
and returns default variable values and indices to those variables."""
function init_variables(variable_model::M) where {T,M <: Model{Varying{T}}}
    deconstructed = var_deconstruct(variable_model)
    initial_p, variable_dxs, p_bounds = init_variables(deconstructed, T)
    return initial_p, variable_dxs, p_bounds
end

function init_variables(deconstructed::Tuple{Type,<:AbstractArray}, T::Type)
    variable_dxs = []
    initial_p = T[]
    p_bounds = Tuple{T,T}[]
    for dx in CartesianIndices(size(deconstructed[2]))
        (typ, val) = deconstructed[2][dx]
        if typ <: Variable
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

function Simulating.time_span(p_search::ParameterSearch)
    time_span(p_search.solver)
end


#region Model Deconstruction & Reconstruction
#########################################################
######### MODEL DECONSTRUCTION & RECONSTRUCTION #########
# Deconstruct a model (which is naturally hierarchical)
# into a flat representation.
# Alter flat representation.
# Reconstruct hierarchical model.
#########################################################


function set_deep_dx!(model_deconstruction::Tuple{Type,Array}, dxs, val)
    tmp = model_deconstruction
    for dx in dxs[1:end-1]
        tmp = tmp[2][dx]
    end
    tmp[2][dxs[end]] = val
end

function deconstruct(v::Variable)
    val = default_value(v)
    return (typeof(val), val)
end

function deconstruct(val::V) where {V <: Union{Number, AbstractString}}
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

function deconstruct(arr::AA) where {AA <: AbstractArray}
    deconstruction = map(deconstruct, arr)
    return (base_type(typeof(arr)), [deconstruction...]) # Remade in case static
end

function deconstruct(tup::Tuple)
    deconstruction = map(deconstruct, tup)
    return (base_type(typeof(tup)), [deconstruction...])
end

function var_deconstruct(val::Variable)
    return (typeof(val), val)
end

function var_deconstruct(val::Union{AbstractString,Number})
    return (typeof(val), val)
end

function var_deconstruct(m::M) where {M <: Parameter}
    deconstruction = Tuple{Type,Any}[]
    for i_field in 1:nfields(m)
        substruct = var_deconstruct(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return (typeof(m), deconstruction)
end

function var_deconstruct(arr::AA) where {AA<:AbstractArray}
    deconstruction = map(var_deconstruct, arr)
    return (base_type(typeof(arr)), [deconstruction...]) # remade incase static
end

function var_deconstruct(tup::Tuple)
    deconstruction = map(var_deconstruct, tup)
    return (base_type(typeof(tup)), [deconstruction...])
end

function base_type(T::Int)
    T
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

function base_type(::Type{Tuple{T,S}}) where {T,S}
    BT = base_type(T)
    BS = base_type(S)
    return Tuple{BT,BS}
end

function base_type(::Type{SArray{TUP,T,N,M}}) where {N,M,T,TUP}
    BT = base_type(T)
    return SArray{TUP,BT,N,M}
end    

function reconstruct(tup::Tuple{Type,<:Union{Number,AbstractString}})
    return tup[2]
end

function reconstruct(tup::Tuple{Type,Array})
    typ, arr = tup
    base_typ = base_type(typ)
    if typ <: Union{AbstractArray, Tuple}
        return base_typ(reconstruct.(arr))
    else
        return base_typ(reconstruct.(arr)...)
    end
end

############################################
##endregion

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
    loss_fn = loss(target, initial_model, solver)
    loss_obj = _build_loss_objective(initial_problem, solver, initial_model.space, loss_fn,
                                                problem_generator)
    result = bboptimize(loss_obj; NumDimensions=length(initial_p),
        MaxSteps=MaxSteps, SearchRange=p_bounds)
    return result
end

function run_search(jl_filename::AbstractString)
    include(jl_filename)
    filecopy(p_search.output, jl_filename, basename(jl_filename))
    return p_search
end

export run_search, model_from_p

end
