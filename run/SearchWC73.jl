module SearchWC73

using Exploration
using WC73
using JLD

function search_WilsonCowan73(jl_filename::String)
	doc"Load jl file and run"
	include(jl_filename)
	p_search::ParameterSearch = load("parameters.jld", "p_search")
	run_search(p_search)
end

end