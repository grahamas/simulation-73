using Modeling, WCM, Meshes, Records, Simulating,
  CalculatedParameters, Analysis, WCMAnalysis,
  WCMConnectivity, WCMNonlinearity, WCMStimulus, WCMTarget
using WCM: WCMSpatial1D
using DifferentialEquations: Euler

if !(@isdefined UV)
  const UV = UnboundedVariable
  const BV = BoundedVariable
  const varying{T} = Union{T,BV{T}}
  const v = Float64
end
T= 4.0
N=1
P=2
simulation = Simulation(
        model = WCMSpatial1D{v,N,P}(;
            pop_names = ["E", "I"],
            α = [1.1, 1.0],
            β = [1.1, 1.1],
            τ = [0.1, 0.18],
            space = Segment{v}(; n_points=101, extent=100),
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
        solver = EulerSolver(;
            T = T,
            dt = 0.01
            ),
        analyses = Analyses{v}(;
          subsampler = SubSampler(;
               space_strides = [4],
               dt = 0.05
               ),
          plot_specs = [
              NeumanTravelingWavePlot(;
                dt = 0.1
              ),
              Animate(;
                fps = 20
                )
              # NonlinearityPlot(;
              #   fn_bounds = (-1,15)),
              # SpaceTimePlot()
              ]
            ),
        output = SingleOutput(;
            root = "/home/grahams/Dropbox/simulation-73/results/",
            simulation_name = "replicate_neuman"
            )
        )

using JLD2

@save "parameters.jld2" simulation
