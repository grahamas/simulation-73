abstract type AbstractOutput end

"""
Output from a single simulation.

    (o::SingleOutput)(write_fn, base_name, args...; kwargs...)

Wraps a `write_fn` to write out to the correct location, "safely."
"""
struct SingleOutput <: AbstractOutput
    root::String
    simulation_name::String
    dir_path::String
    safe_writer::Function
end

function SingleOutput(; root=nothing, simulation_name=nothing, dir_prefix=nothing)
    dir_path = directory(root, simulation_name, dir_prefix)
    writer = make_writer(dir_path)
    SingleOutput(root, simulation_name, dir_path, writer)
end

function (o::SingleOutput)(write_fn::Function, base_name::String, args...; kwargs...)
    o.safe_writer(write_fn, base_name, args...; kwargs...)
end


raw"""
    directory(root, simulation_name, prefix)

Construct a directory name in a directory `simulation_name` at `root`, whose name is `$(prefix)_$(now)` where `$(now)` is the current timestamp.
"""
function directory(root, simulation_name, prefix)
    nowstr = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS.s")
    basename = "$(prefix)_$(nowstr)"
    dir_name = joinpath(root, simulation_name, nowstr)
    mkpath(dir_name)
    return dir_name
end

# FIXME: Reason for prefix unclear
"""
    full_path(dir_name, base_name, prefix="")

Construct the full path to a file, given the containing directory `dir_name`, basename `base_name`, and an optional prefix for the basename.
"""
function full_path(dir_name::String, base_name::String, prefix::String="")
    prefixed_name = join([x for x in [prefix, base_name] if length(x) > 0], "_")
    full_path = joinpath(dir_name, prefixed_name)
    return full_path
end

"""
    make_writer(dir_path)

Return a functor to "safely" wrap a file-writing function.

Note "safe" here means it won't clobber an existing file.
"""
function make_writer(dir_path::String)
    function safe_write_fn(write_fn::Function, base_name, args...; kwargs...)
        fp = full_path(dir_path, base_name)
        if !(isfile(fp))
            @info "Writing $fp"
            write_fn(fp, args...; kwargs...)
        else
            @warn("Tried to write existing file: $fp")
        end
    end
    return safe_write_fn
end

"""
    filecopy(output, source_path, dest_base_name)

Copy a fully specified source `source_path` to `dest_base_name` via `output`.
"""
function filecopy(output::AbstractOutput, source_path, dest_base_name)
    cp(source_path, joinpath(output.dir_path, dest_base_name))
end
