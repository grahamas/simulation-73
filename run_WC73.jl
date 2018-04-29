if !("./run" in LOAD_PATH)
    push!(LOAD_PATH, "./run")
    import RunWC73
else
    reload("RunWC73.jl")
end

changes = ["heightened_EI", "original"]
nonlinearity_exts = ["", "_sigmoid_diff"]#, "_sech2"]
ext = ".json"

param_files = [joinpath("params", change .* nonlinearity .* ext) for change in changes, nonlinearity in nonlinearity_exts][:]

RunWC73.run_WilsonCowan73_trial(param_files, Dict(:model =>
                                                  Dict(
                                                      #o:τ => RowVector([0.18, 0.1]),
                                                      #:connectivity => Dict(:spreads => [3.0 2.7; 2.7 2.5]),
                                                      :stimulus => Dict(
                                                          :name => "smooth_bump",
                                                          :strength => 10.0,
                                                          :steepness => 1.0
                                                      ),
                                                      :nonlinearity => Dict(
                                                          :args => Dict(
                                                              :a => RowVector([1.3, 0.8]),
                                                              :θ => RowVector([8.0, 2.6]))))))
