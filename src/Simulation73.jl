module Simulation73

using Markdown # for doc_str
using CalculatedTypes
import CalculatedTypes: Calculated, update
using DifferentialEquations, DiffEqParamEstim
using BlackBoxOptim, Optim
using StaticArrays
import BlackBoxOptim
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve, Euler, ODEProblem
using Dates
using RecipesBase
using Parameters
using Plots

# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"
# using Plots

# "variables.jl"
export AbstractVariable, UnboundedVariable, BoundedVariable,
	default_value, bounds, pops, MaybeVariable,
	AbstractParameter

# "records.jl"
export SingleOutput, ExperimentOutput, filecopy, AbstractOutput

# "meshes.jl"
export DistanceMatrix, CalculatedDistanceMatrix, get_origin,
	AbstractSpace, Pops, AbstractLine,
	Segment, Pops, Circle

# "subsampling.jl" (note: should probably be meshed with meshes)
export scalar_to_idx_window, subsampling_Δidx, subsampling_idxs,
	subsampling_time_idxs, subsampling_space_idxs, Subsampler

# "analysing.jl"
export AbstractPlotSpecification, AbstractSpaceTimePlotSpecification, Analyses,
	output_name, plot_and_save, analyse, subsample, subsampling_idxs

# "targets.jl"
export AbstractTarget, target_loss

# "simulating.jl"
export Model, initial_value, space_arr,
	run_simulation, Simulation, write_params, Solver,
    time_span, time_arr, save_dt, save_dx, save_idxs, generate_problem

# "exploring.jl"
export ParameterSearch, make_problem_generator, update_from_p!,
	make_calculated_function

export run_search, model_from_p, result_simulation

include("variables.jl")
include("records.jl")
include("destructuring.jl")
include("meshes.jl")
include("subsampling.jl")
include("simulating.jl")
include("targets.jl")
include("exploring.jl")
include("analysing.jl")

end
