module Records

using Parameters
using JLD2
using Memoize
using Dates
using Logging

abstract type AbstractOutput end

mutable struct SingleOutput <: AbstractOutput
    root::String
    simulation_name::String
    dir_path::String
    safe_writer::Function
end
function SingleOutput(; root=nothing, simulation_name=nothing)
    dir_path = directory(root, simulation_name)
    writer = make_writer(dir_path)
    SingleOutput(root, simulation_name, dir_path, writer)
end

# @with_kw struct ExperimentOutput <: AbstractOutput
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

function directory(root, simulation_name)
    nowstr = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS.s")
    dir_name = joinpath(root, simulation_name, nowstr)
    mkpath(dir_name)
    return dir_name
end

function full_path(dir_name::String, base_name::String, prefix::String="")
    prefixed_name = join([x for x in [prefix, base_name] if length(x) > 0], "_")
    full_path = joinpath(dir_name, prefixed_name)
    return full_path
end

function make_writer(dir_path::String)
    function safe_write_fn(write_fn::Function, base_name, args...; kwargs...)
        fp = full_path(dir_path, base_name)
        if !(isfile(fp))
            @info "Writing $fp"
            write_fn(fp, args...; kwargs...)
        else
            warn("Tried to write existing file: $fp")
        end
    end
    return safe_write_fn
end

function filecopy(output::AbstractOutput, source_path, dest_base_name)
    cp(source_path, joinpath(output.dir_path, dest_base_name))
end

function (o::SingleOutput)(write_fn::Function, base_name::AbstractString, args...; kwargs...)
    o.safe_writer(write_fn, base_name, args...; kwargs...)
end

function (o::SingleOutput)(obj; base_name::AbstractString, write_fn)
    return nothing
end
function required_modules()
    error("undefined.")
end

export AbstractOutput, SingleOutput, ExperimentOutput, filecopy

end