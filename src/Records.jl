module Records

using Parameters
using JLD2
using Memoize

abstract type Output end

@with_kw struct SingleOutput <: Output
    root::String
    simulation_name::String
end


# @with_kw struct ExperimentOutput <: Output
#     root::String
#     simulation_name::String
#     mod_name::String
#     experiment_name::String
# end

# * File ops

# function make_experiment_output_folder(root, simulation_name, mod_name, experiment_name)
#     @assert length(mod_name) > 0
#     #now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
#     #experiment_dir_name = join([experiment_name, now], "_")
#     dir_name = joinpath(root, experiment_name)
#     mkpath(dir_name)
#     return (dir_name, if (length(simulation_name) > 0) join([simulation_name, mod_name], "_") else mod_name end)
# end

@memoize function directory(output::SingleOutput)
    now = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.s")
    dir_name = joinpath(output.root, output.simulation_name, now)
    mkpath(dir_name)
    return dir_name
end

function make_writer(dir_name::String, prefix::String="")
    function safe_write_fn(write_fn::Function, base_name)
        prefixed_name = join([x for x in [prefix, base_name] if length(x) > 0], "_")
        full_path = joinpath(dir_name, prefixed_name)
            if !(isfile(full_path))
                write_fn(full_path)
            else
                warn("Tried to write existing file: $full_path")
            end
    end
    return safe_writer
end

@memoize function make_writer(output::SingleOutput)
    root = output.root
    simulation_name = output.simulation_name
    dir_name = directory(output)
    make_writer(dir_name)
end

(o::SingleOutput)(write_fn::Function, base_name::AbstractString) = begin
    output_fn = make_writer(o)
    output_fn(write_fn, base_name)
end

# function make_write_fn(output::ExperimentOutput)
#     root = output.root
#     simulation_name = output.simulation_name
#     mod_name = output.mod_name
#     experiment_name = output.experiment_name
#     dir_name, prefix = make_experiment_output_folder(root, simulation_name, mod_name, experiment_name)
#     make_write_fn(dir_name, prefix)
# end

function required_modules()
    error("undefined.")
end

function write_object(output::Output, file_name::String, object_name::String, object)
    object_write(path) = jldopen(path, "w") do file
        #addrequire.(file, required_modules(typeof(object)))
        file[object_name] = object
    end
    write_fn(output)(object_write, file_name)
end

export Output, SingleOutput, ExperimentOutput
export write_object, required_modules, write_fn

end