module Records

using Parameters
using JLD

abstract type Output end

@with_kw struct SingleSimulation <: Output
    root::String
    simulation_name::String
end

# * File ops
function make_individual_output_folder(root, simulation_name, mod_name)
    now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
    if length(mod_name) > 0
        now = join([mod_name, now], "_")
    end
    dir_name = joinpath(root, simulation_name, now)
    mkpath(dir_name)
    return dir_name
end

function make_experiment_output_folder(root, simulation_name, mod_name, experiment_name)
    @assert length(mod_name) > 0
    #now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
    #experiment_dir_name = join([experiment_name, now], "_")
    dir_name = joinpath(root, experiment_name)
    mkpath(dir_name)
    return (dir_name, if (length(simulation_name) > 0) join([simulation_name, mod_name], "_") else mod_name end)
end


function make_write_fn(dir_name::String, prefix::String="")
     function safe_write_fn(write_fn, base_name)
        prefixed_name = join([prefix, base_name], "_")
        full_path = joinpath(dir_name, prefixed_name)
        if !(isfile(full_path))
   	     write_fn(full_path)
	else
	     warn("Tried to write existing file: $full_path")
	end
    end
    return safe_write_fn
end
e
function make_write_fn(output::SingleOutput)
    root = output.root
    simulation_name = output.simulation_name
    dir_name = make_individual_output_folder(root, simulation_name)
    make_write_fn(dir_name)
end

function make_write_fn(output::ExperimentOutput)
    root = output.root
    simulation_name = output.simulation_name
    mod_name = output.mod_name
    experiment_name = output.experiment_name
    dir_name, prefix = make_experiment_output_folder(root, simulation_name, mod_name, experiment_name)
    make_write_fn(dir_name, prefix)
end

function write_params(safe_write_fn; params...)
    base_name = "parameters.json"
    safe_write_fn((path) -> write(path, JSON.json(params,4)), base_name)
end

function write_params(safe_write_fn, sim::Simulation)
    base_name = "parameters.jld"
    safe_write_fn(base_name) do path
        jldopen(path, "w") do file
            addrequire(file, Simulation)
            write(file, "sim", sim)
        end
    end
end
