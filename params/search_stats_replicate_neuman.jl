using Modeling, WCM, Meshes, Records, Simulating,
  CalculatedParameters, Analysis, WCMAnalysis,
  WCMConnectivity, WCMNonlinearity, WCMStimulus,
  WCMTarget
using WCM: WCMSpatial1D
using DifferentialEquations: Euler

const v = Varying{Float64}

stop_time = 3.0
T= 4.0
N=1
P=2
@show v
p_search = ParameterSearch(;
        varying_model = WCMSpatial1D{v,N,P}(;
            pop_names = ["E", "I"],
            α = v[1.1, 1.0],
            β = v[1.1, 1.1],
            τ = v[BV(0.1, (0.05,0.25)), 0.18],
            space = PopSegment{v,P}(; n_points=301, extent=100.0),
            nonlinearity = pops(SigmoidNonlinearity{v}; a=v[1.2, 1.0],
                                                        θ=v[2.6, 8.0]),
            stimulus = pops(NoisySharpBumpStimulus{v}; strength=v[1.2, 1.2],
                                                   window=Tuple{v,v}[(0.0,0.55), (0.0,0.55)],
                                                   width=v[2.81, 2.81],
                                                   SNR=v[80.0, 80.0]),
            connectivity = pops(ShollConnectivity{v};
                amplitude = v[16.0 -18.2;
                             27.0 -4.0],
                spread = v[2.5 2.7;
                           2.7 2.5])
            ),
        solver = Solver(;
            stop_time = 2.0,
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
            simulation_name = "plotting_tests"
            ),
        target=MatchExample(
            file_name="data/wave_example.jld2"
          ),
        MaxSteps=100
        )

using JLD2
@save "parameters.jld2" p_search
