module SearchWC73

using Exploration
using WC73
using JLD2

"Load jl file and search parameter space"
function search_WilsonCowan73(jl_filename::String)
	include(jl_filename)
	@load "parameters.jld2" p_search
	result_sim, result_soln = run_search(p_search)
	return (p_search, result_sim, result_soln)
end

end