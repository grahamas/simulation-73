module RunWC73

using WC73

first(pg::ParamGenerator) = next(pg, start(pg))[1]
function run_WilsonCowan73_trial(json_filename::String, varying_mods::Array,
                                 constant_mods=nothing::Union{Dict, Void})
    base_params = WC73.load_WilsonCowan73_parameters(json_filename, constant_mods)
    #run_WilsonCowan73_trial(first(ParamGenerator(base_params, varying_mods)))
    pmap(run_WilsonCowan73_trial, ParamGenerator(base_params, varying_mods))
end

end
