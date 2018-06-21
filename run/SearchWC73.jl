module SearchWC73

using Exploration
using WC73

function search_WilsonCowan73(jl_filename::String)
	doc"Load jl file and run"
	include(jl_filename) # Defines function "make_parameters"
	p_search::ParameterSearch = make_parameters()
	run(p_search)
end

end