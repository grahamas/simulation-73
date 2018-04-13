push!(LOAD_PATH, ".")

using WC73: load_WilsonCowan73_parameters, solve_WilsonCowan73
using WC73Analysis: analyse_WilsonCowan73_solution

function run_WilsonCowan73_trial(json_filename::String, modifications=nothing::Union{Dict, Void})
    all_params = load_WilsonCowan73_parameters(json_filename, modifications)
    solution = solve_WilsonCowan73(; all_params...)
    analyse_WilsonCowan73_solution(solution; all_params...)
end
