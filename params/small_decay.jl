const UV = UnboundVariable

function make_target(model::WC73)
    segment = calculate(model.space).value
    c = 10
    wave_fn(t) = @. exp(-2.5t) * sech(segment - c * t)
    return wave_fn.(1:3)
end


function make_parameters()
    ParameterSearch(
        model = WC73(
            α = [1.1, 1.0],
            β = [1.1, 1.1],
            τ = [0.1, 0.18],
            space = Segment(n_points=401, extent=100.5)
            nonlinearity = pops(SigmoidNonlinearity; a=[1.2, 1.0], θ=[2.6, 8.0])
            stimulus = pops(SharpBumpStimulus; strength=[2.4,2.4], duration=[1.0,1.0], width=[3.0,3.0])
            connectivity = pops(ShollConnectivity;
                                amplitude = [UV(16.0) -18.2;
                                             UV(27.0) -4.0],
                                spread = [2.5 2.7;
                                          2.7 2.5])
        ),
        solver = Solver(
            T = 6,
            params = Dict(
                :dt => 0.001,
                #:alg_hints => [:stiff]
            )
        )
        analyses =  Dict(
	    "pop_names" => ["E", "I"],
	    "down_sampling" => Dict(
	        "spatial_stride" => 2,
	        "temporal_stride" => 20
	    ),
	    "activity_gif" => Dict(
	        "file_name" => "activity.gif",
	        "disable" => 0,
	        "fps" => 20
	    )
        ),
        output = Output(
            root = "/home/grahams/Dropbox/Research/simulation-73/results/",
            simulation_name = ""
        ),
    target = Target(
        factory=make_target
    )
)
end
