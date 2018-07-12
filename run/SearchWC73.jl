module SearchWC73

macro print_using(mod)
	quote
        #println("Using ", $(string(mod)))
        using $mod
        #println("... done using ", $(string(mod)))
    end
end
@print_using Exploration
@print_using WC73
using JLD

function search_WilsonCowan73(jl_filename::String)
	doc"Load jl file and run"
	include(jl_filename)
	p_search::ParameterSearch = load("parameters.jld", "p_search")
	run_search(p_search)
end

end