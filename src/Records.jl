module Records

using Parameters
using JLD2
using Memoize
using Dates

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
#     #nowstr = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS.s")
#     #experiment_dir_name = join([experiment_name, nowstr], "_")
#     dir_name = joinpath(root, experiment_name)
#     mkpath(dir_name)
#     return (dir_name, if (length(simulation_name) > 0) join([simulation_name, mod_name], "_") else mod_name end)
# end

@memoize function directory(output::SingleOutput)
    nowstr = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS.s")
    dir_name = joinpath(output.root, output.simulation_name, nowstr)
    mkpath(dir_name)
    return dir_name
end

function full_path(dir_name::String, base_name::String, prefix::String="")
    prefixed_name = join([x for x in [prefix, base_name] if length(x) > 0], "_")
    full_path = joinpath(dir_name, prefixed_name)
    return full_path
end

function safe_write_fn(write_fn::Function, base_name, args...; kwargs...)
    if !(isfile(full_path))
        write_fn(full_path, args...; kwargs...)
    else
        warn("Tried to write existing file: $full_path")
    end
end

@memoize function make_writer(output::SingleOutput)
    root = output.root
    simulation_name = output.simulation_name
    dir_name = directory(output)
    make_writer(dir_name)
end

function (o::SingleOutput)(write_fn::Function, base_name::AbstractString, args...; kwargs...)
    output_fn = make_writer(o)
    output_fn(write_fn, base_name, args...; kwargs...)
end

function (o::SingleOutput)(obj; base_name::AbstractString, write_fn)

function required_modules()
    error("undefined.")
end

export Output, SingleOutput, ExperimentOutput

end