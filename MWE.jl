
using DifferentialEquations, DiffEqParamEstim
using BlackBoxOptim

function sigmoid_fn(x, a, theta)
    return @. (1.0 / (1 + exp(-a * (x - theta))))
end

function make_ode_fn(the_p::Array{T,1}) where T
	τ = the_p[1]
	function ode_fn!(dA::Array{T,1}, A::Array{T,1}, p::Union{Array{T,1},Nothing}, t::T)
		@. dA = (-0.1 * A + 1.1 * (1.0 - A) * sigmoid_fn(3.0 * A + 0.5, 1.0, 8.0)) / τ
	end
	return ode_fn!
end

function problem_generator(prob, new_p::Array{T,1}) where T
	fn = make_ode_fn(new_p)
	return ODEProblem(fn, [0.0; 0.1], (0.0, 1.0), new_p)
end

function run_search()
	target = zeros(Float64,2,11)
	initial_problem = problem_generator(nothing, [1.0])
	loss_fn = L2Loss(0.0:0.1:1.0, target)
	loss_obj = build_loss_objective(initial_problem, Euler(), loss_fn;
		prob_generator=problem_generator, dt=0.1)
	result = bboptimize(loss_obj; NumDimensions=1,
		MaxSteps=10, SearchRange=[(0.1,4.0)])
	return result
end

function solve_from_result(result)
	p = best_candidate(result)
	fn = make_ode_fn(p)
	prob = ODEProblem(fn, [0.1, 0.0], (0.0, 1.0), nothing)
	solve(prob, Euler(), dt=0.1)
end