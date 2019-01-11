module Records

using Parameters
using JLD2
using Memoize
using Dates
using Logging

abstract type AbstractOutput end

@with_kw mutable struct SingleOutput <:AbstractOutput
    root::String
    simulation_name::String
end

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

function make_writer(dir_name::AbstractString)
    function safe_write_fn(write_fn::Function, base_name, args...; kwargs...)
        fp = full_path(dir_name, base_name)
        if !(isfile(fp))
            @info "Writing $fp"
            write_fn(fp, args...; kwargs...)
        else
            warn("Tried to write existing file: $fp")
        end
    end
    return safe_write_fn
end

@memoize function make_writer(output::SingleOutput)
    dir_name = directory(output)
    make_writer(dir_name)
end

function (o::SingleOutput)(write_fn::Function, base_name::AbstractString, args...; kwargs...)
    output_fn = make_writer(o)
    output_fn(write_fn, base_name, args...; kwargs...)
end

function (o::SingleOutput)(obj; base_name::AbstractString, write_fn)
    return nothing
end
function required_modules()
    error("undefined.")
end

export AbstractOutput, SingleOutput, ExperimentOutput

end