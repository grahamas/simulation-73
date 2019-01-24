using Modeling, WCM, Meshes, Records, Simulating,
  CalculatedParameters, Analysis, WCMAnalysis,
  WCMConnectivity, WCMNonlinearity, WCMStimulus
using WCM: WCMSpatial1D
using DifferentialEquations: Euler

if !(@isdefined UV)
  const UV = UnboundedVariable
  const BV = BoundedVariable
  const varying{T} = Union{T,BV{T}}
  const v = Float64
end
stop_time = 3.0
T= 4.0
N=1
P=2
simulation = Simulation(
        model = WCMSpatial1D{v,N,P}(;
            pop_names = ["E", "I"],
            α = v[1.1, 1.0],
            β = v[1.1, 1.1],
            τ = v[0.1, 0.18],
            space = PopSegment{v,P}(; n_points=301, extent=100.0),
            nonlinearity = pops(SigmoidNonlinearity{v}; a=[1.2, 1.0],
                                                        θ=[2.6, 8.0]),
            stimulus = pops(NoisySharpBumpStimulus{v}; strength=[1.2, 1.2],
                                                   window=[(0.0,0.55), (0.0,0.55)],
                                                   width=[2.81, 2.81],
                                                   SNR=[80.0, 80.0]),
            connectivity = pops(ShollConnectivity{v};
                amplitude = [16.0 -18.2;
                             27.0 -4.0],
                spread = v[2.5 2.7;
                           2.7 2.5])
            ),
        solver = Solver(;
            stop_time = 3.0,
            dt = 0.01,
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
              PeakTravelingWavePlot(;
                dt = 0.1
              )
              ]
            ),
        output = SingleOutput(;
            root = "/home/grahams/Dropbox/simulation-73/results/",
            simulation_name = "replicate_neuman"
            )
        )

using JLD2

@save "parameters.jld2" simulation
