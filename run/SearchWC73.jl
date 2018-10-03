module SearchWC73

using Exploration
using WC73
using FileIO

"Load jl file and run"
function search_WilsonCowan73(jl_filename::String)
	include(jl_filename)
	p_search = load("parameters.jld2", "p_search")
	run_search(p_search)
end

end