module SearchWC73

using Exploration
using WC73
using JLD2

function search_WilsonCowan73(jl_filename::String)
	doc"Load jl file and run"
	include(jl_filename)
	p_search::ParameterSearch = load("parameters.jld2", "p_search")
	result_sim, result_soln = run_search(p_search)
	return (p_search, result_sim, result_soln)
end

end