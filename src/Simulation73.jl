module Simulation73

using DrWatson

using Markdown # for doc_str
using DifferentialEquations, DiffEqParamEstim
using BlackBoxOptim, Optim
using StaticArrays
import BlackBoxOptim
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve, Euler, ODEProblem
using RecipesBase
using Parameters
using Plots

# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"
# using Plots

abstract type AbstractParameter{T} end
DrWatson.default_allowed(c::AbstractParameter) = (Real, String, Symbol, AbstractParameter)

# "variables.jl"
export AbstractVariable, UnboundedVariable, BoundedVariable,
	default_value, bounds, pops, MaybeVariable,
	AbstractParameter

# "meshes.jl"
export Pops, Segment, Circle,
	 Grid, Torus

export AbstractSpace, DistanceMatrix, CalculatedDistanceMatrix, get_origin,
	 get_distances, one_pop, one_pop_zero, coordinates

# "subsampling.jl" (note: should probably be meshed with meshes)
export scalar_to_idx_window, subsampling_Δidx, subsampling_idxs,
	subsampling_time_idxs, subsampling_space_idxs, Subsampler

# "analysing.jl"
export AbstractPlotSpecification, AbstractSpaceTimePlotSpecification, Analyses,
	output_name, plot_and_save, analyse, subsample, subsampling_idxs

# "targets.jl"
export AbstractTarget, target_loss

export execute

# "simulating.jl"
export AbstractModel, Solver, Simulation, Execution,
	initial_value, time_span, save_dt, save_dx,
	generate_problem, solve, run_simulation,
	make_mutators, make_system_mutator

# "exploring.jl"
export Search, SearchExecution, make_problem_generator, search, run_search

include("helpers.jl")
include("deconstructing.jl")
include("variables.jl")
include("records.jl")
include("meshes.jl")
include("subsampling.jl")
include("simulating.jl")
include("targets.jl")
include("exploring.jl")
include("analysing.jl")

end
