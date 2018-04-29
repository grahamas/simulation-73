module RunWC73

reload("WC73")
using WC73

push!(LOAD_PATH, ".")
reload("Solitons")
using Solitons

using Plots; pyplot()
using Optim

# * Run simulation
function run_WilsonCowan73_trial(json_filename::String, modifications=nothing::Union{Dict, Void})
    all_params = WC73.load_WilsonCowan73_parameters(json_filename, modifications)
    solution = WC73.solve_WilsonCowan73(; all_params...)
    WC73.analyse_WilsonCowan73_solution(solution; all_params...)
end

function run_WilsonCowan73_trial(l_json_filename::Array{String,1}, modifications=nothing::Union{Dict,Void})
    @sync @parallel for json_filename in l_json_filename
        run_WilsonCowan73_trial(json_filename, modifications)
    end
end


# * Example solitons

# ** make_example_solitons

function make_example_solitons(T::DataType)
    mesh = SpaceMesh(T,[Dict(:extent=>100,:N=>281)])
    time_len = parse(T, "15")
    dt = parse(T, "0.01")
    make_example_solitons(T, mesh, time_len, dt)
end

function make_example_solitons(ValueT::DataType,
                               mesh::SpaceMesh{DistT},
                               time_len::TimeT,
                               dt::TimeT) where {DistT<:Real, TimeT<:Real}
    function make_example_soliton(soliton_fn::Function)
        example = Dict{Symbol, Any}()
        example[:soliton] =
            apply_through_time(soliton_fn, mesh, time_len, dt)
        example[:mesh] = mesh
        example[:time_len] = time_len
        example[:dt] = dt
        return example
    end

    parse_param(param_str) = parse(ValueT, param_str)

    example_solitons = Dict()

    simple_velocity = parse_param("20")
    simple_x0 = parse_param("0")
    simple_soliton_fn = make_simple_soliton(simple_velocity, simple_x0)
    example_solitons[:simple] = make_example_soliton(simple_soliton_fn)

    td_param_strs = ["10", "20", "0.5", "0.5", "1", "0"]
    td_params = parse_param.(td_param_strs)
    td_soliton_fn = make_time_dependent_soliton(td_params...)
    example_solitons[:td] = make_example_soliton(td_soliton_fn)

    function mesh_percentile(mesh, p)
        dim = mesh.dims[1]
        dim[floor(Int,length(dim)*p)]
    end
    left_velocity = parse_param("20")
    right_velocity = parse_param("-20")
    left_soliton_fn = make_simple_soliton(left_velocity, mesh_percentile(mesh,0.25))
    right_soliton_fn = make_simple_soliton(right_velocity,  mesh_percentile(mesh,0.75))
    simple_collision_fn(x, t) = left_soliton_fn(x,t) + right_soliton_fn(x,t)
    example_solitons[:simple_collision] = make_example_soliton(simple_collision_fn)

    return example_solitons
end

# ** plot_soliton and plot_example_soliton

function plot_soliton(save_path; soliton=nothing, mesh=nothing, time_len=nothing, dt=nothing)
    @assert all([soliton, mesh, time_len, dt] .!= nothing)
    plot_solution_surface(soliton, mesh, time_len, dt,
                          save=save_path)
end


function plot_example_solitons(examples, save_path)
    for (name, soliton_dct) in examples
        plot_soliton(joinpath(save_path, (String(name) * "_soliton.png")); soliton_dct...)
    end
end

# * Local peak detection

function amplitude_width_offset(values, mesh, max_dx, half_window::Int=1)
    @assert length(mesh.dims) == 1
    around_max = values[max_dx-half_window:max_dx+half_window]
    locs = mesh.dims[1][max_dx-half_window:max_dx+half_window]

    max_guess = around_max[half_window+1]
    upper_bound_max = max_guess + maximum(abs.(
        max_guess .- around_max[half_window:half_window+2][[1,3]]))
    guess = [max_guess, 1, locs[half_window+1]]
    lower = [max_guess - 10*eps(), 0, locs[half_window]]
    upper = [upper_bound_max, Inf, locs[half_window+2]]
    return around_max, locs, guess, lower, upper
end

function get_peak_maximum_parabola_fit(args...)
    parabola(x,p) = -p[3] .* ((x .- p[2]) .^ 2) .+ p[1]
    space_args = amplitude_width_offset(args...)
    return get_peak_maximum_fit(parabola, space_args...)
end

function get_peak_maximum_sech2_fit(args...)
    sech2(x,p) = p[1] .* (sech.(p[2] .* (x .- p[3])) .^ 2)
    jac_sech2(x,p) = [sech.(p[2] .* (x .- p[3])) .^ 2,
                      -2 .* sech2(x,p) .* tanh.(p[2] .* (x - p[3])) .* (x - p[3]),
                      2 .* p[2] .* sech2(x,p) .* tanh.(p[2] .* (x - p[3])) .* (x - p[3])] # Verify dimensions
    space_args = amplitude_width_offset(args...)
    return get_peak_maximum_fit(sech2, space_args...)
end

function get_peak_maximum_fit(fn, around_max, locs, guess, lower, upper)
    # fn(x,p)
    err_fn(p) = norm(fn(locs,p) - around_max)

    fit = optimize(err_fn, guess, lower, upper, Fminbox{GradientDescent}())

    amplitude = Optim.minimizer(fit)[1]
    if amplitude > 10.5
        plot(locs, around_max, seriestype=:scatter)
        xs = linspace(locs[1], locs[end], 100)
        plot!(xs, fn(xs, Optim.minimizer(fit)))
        savefig("/home/grahams/Dropbox/target/example_fit_$(amplitude).png")
    end
    return amplitude
end


function get_peak_maximum(values, mesh, max_dx, args...)
    values[max_dx]
end

function plot_peak_magnitude{ValueT,TimeT}(save_path, fit_half_window=1;
                                           soliton::Array{ValueT}=nothing,
                                           mesh=nothing, time_len::TimeT=nothing,
                                           dt::TimeT=nothing)
    @assert all([soliton, mesh, time_len, dt] .!= nothing)
    @assert eltype(soliton) == ValueT
    timepoints = 0:dt:time_len
    @assert eltype(timepoints) == TimeT
    maxima = zeros(ValueT, size(timepoints))
    plot_soliton((save_path * "_soliton.png"); soliton=soliton, mesh=mesh, time_len=time_len, dt=dt)
    for (t_dx, t) in enumerate(timepoints)
        # try
            t_maxima_dx = find_local_maxima(soliton[:,t_dx])
            maxima[t_dx] = get_peak_maximum_sech2_fit(soliton[:,t_dx], mesh, t_maxima_dx[1], fit_half_window)
        # catch err
        #     println(err)
        #     maxima[t_dx] = 0
        # end
    end
    println(maxima)
    plot(timepoints, maxima)
    savefig(save_path * "_peak_maxima.png")
end


function plot_peak_magnitudes(examples, save_path, fit_half_window=1)
    for (name_sym, soliton_dct) in examples
        plot_peak_magnitude(joinpath(save_path, String(name_sym)), fit_half_window; soliton_dct...)
    end
end

end
