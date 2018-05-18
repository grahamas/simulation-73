if !("./run" in LOAD_PATH)
    push!(LOAD_PATH, "./run")
    import RunWC73
else
    reload("RunWC73.jl")
end

# changes = ["heightened_EI", "original"]
# nonlinearity_exts = ["_sigmoid_diff"]#, "_sech2"]
# ext = ".json"

# param_files = [joinpath("params", change .* nonlinearity .* ext) for change in changes, nonlinearity in nonlinearity_exts][:]

param_file = "params/original.json"

RunWC73.run_WilsonCowan73_trial(param_file, [
                                #  ((:model, :β, 2), [0.8, 1.1, 1.5, 2.0, 3.0], "βI"),
                                #   ((:model, :connectivity, :amplitudes, (1,1)), 2.0:1.0:10.0, "EE"),
                                #   ((:model, :connectivity, :amplitudes, (2,1)), 1.5:1.0:10.5, "EI"),
      #                           #   ((:model, :connectivity, :amplitudes, (2,2)), [-4.0, -2.0], "II"),
     ((:model, :stimulus, :strength), 0.3647807:0.00000001:0.3647808, "strength")
     #((:model, :stimulus, :name), ["sharp_bump", "two_sharp_bumps"], "stim"),
     #((:model, :stimulus, :duration), 0.1:0.2:0.3, "duration")
    # #  ((:model, :nonlinearity, :args, :width, 1), [3.0, 5.0], "widthE")],
    ],
                                Dict(:output => Dict(:experiment_name => "poster_strength_finerer8"),
                                    :model =>
                                     Dict(
                                         #     :τ => RowVector([0.18, 0.1]),
                                         #     #:connectivity => Dict(:spreads => [3.0 2.7; 2.7 2.5]),
                                         #     # :stimulus => Dict(
                                         #     #     #:name => "sharp_bump",
                                         #     #     #:strength => 3.0,
                                         #     #     :duration => 2.0
                                         #     #     #:steepness => 1.0
                                             ),
                                     :analyses => Dict(
                                          :activity_gif => Dict(
                                              :disable => 0
                                         ),
                                         :heatmap => Dict(
                                             :name => "heatmap.png"
                                         ),
                                         :nonlinearity => Dict(
                                             :disable => 1
                                         )#,

                                     )
                                     )
                                )
