module Simulation73

using DifferentialEquations, DiffEqParamEstim
using BlackBoxOptim
using StaticArrays
import BlackBoxOptim: OptimizationResults
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve, Euler, ODEProblem

# "meshes.jl"
export DistanceMatrix, CalculatedDistanceMatrix,
	Space, PopSegment, PopSpace, PopSegment

# "subsampling.jl" (note: should probably be meshed with meshes)
export scalar_to_idx_window, subsampling_Δidx, subsampling_idxs 

# "modeling.jl"
export Model, initial_value, space_arr

# "simulating.jl"
export run_simulation, Simulation, write_params, Solver,
    time_span, time_arr, save_dt, save_dx, save_idxs, generate_problem

# "analysing.jl"
export AbstractPlotSpecification, AbstractSpaceTimePlotSpecification, Analyses, 
	output_name, plot_and_save, analyse

# "targets.jl"
export Target, TargetFactory, LossFunction, loss

# "exploring.jl"
export ParameterSearch, Varying, make_problem_generator


include("meshes.jl")
include("subsampling.jl")
include("modeling.jl")
include("simulating.jl")
include("analysing.jl")
include("targets.jl")
include("exploring.jl")

end