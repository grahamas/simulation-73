using Modeling, WCM, Meshes, Records, Simulating,
  CalculatedParameters, Analysis, WCMAnalysis,
  WCMConnectivity, WCMNonlinearity, WCMStimulus,
  WCMTarget
using WCM: WCMSpatial1D
using DifferentialEquations: Euler

const v = Varying{Float64}

N=1
P=2
@show v
p_search = ParameterSearch(;
        varying_model = WCMSpatial1D{v,N,P}(;
            pop_names = ["E", "I"],
            α = v[1.1, 1.0],
            β = v[1.1, 1.1],
            τ = v[0.18, 0.1], #REVERSED
            space = PopSegment{v,P}(; n_points=301, extent=100.0),
            nonlinearity = pops(SigmoidNonlinearity{v}; a=v[BV(1.2, (0.9, 1.3)), BV(1.0, (0.9, 1.3))],
                                                        θ=v[8.0, 2.6]),
            stimulus = pops(NoisySharpBumpStimulus{v}; strength=v[1.2, 1.2],
                                                   window=Tuple{v,v}[(0.0,0.55), (0.0,0.55)],
                                                   width=v[2.81, 2.81],
                                                   SNR=v[80.0, 80.0]),
            connectivity = pops(ShollConnectivity{v};
                amplitude = v[BV(16.0, (8.0,24.0)) BV(-18.2, (-30.0,-15.0));
                             BV(27.0, (15.0,30.0)) BV(-4.0,(-12.0,-2.0))],
                spread = v[2.5 2.7;
                           2.7 2.5])
            ),
        solver = Solver(;
            stop_time = 1.8,
            dt = 0.005,
            space_save_every=1,
            time_save_every=1,
            algorithm=Euler()
            ),
        analyses = Analyses(;
          plot_specs = [
              # Animate(;
              #   fps = 20
              #   ),
              # NonlinearityPlot(;
              #   fn_bounds = (-1,15)
              #   ),
              # SpaceTimePlot(),
              SubsampledPlot(
                plot_type=WaveStatsPlot,
                time_subsampling=Dict(
                    :Δsubsampled => 0.01,
                    :scalar_window => (1.2, 1.8)
                  ),
                space_subsampling=Dict(
                    :scalar_window => (5.0,Inf)
                  )
              )
              ]
            ),
        output = SingleOutput(;
            root = "/home/grahams/Dropbox/simulation-73/results/",
            simulation_name = "broad_search/reversing_e_i"
            ),
        target=MatchExample(
            file_name="data/wave_example.jld2"
          ),
        MaxSteps=10000
        )

using JLD2
@save "parameters.jld2" p_search
