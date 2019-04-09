abstract type AbstractOutput end

"Output from a single simulation."
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

function directory(root, simulation_name, prefix)
    nowstr = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS.s")
    basename = "$(prefix)_$(nowstr)"
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
            @warn("Tried to write existing file: $fp")
        end
    end
    return safe_write_fn
end

function filecopy(output::AbstractOutput, source_path, dest_base_name)
    cp(source_path, joinpath(output.dir_path, dest_base_name))
end

function (o::SingleOutput)(write_fn::Function, base_name::String, args...; kwargs...)
    o.safe_writer(write_fn, base_name, args...; kwargs...)
end
