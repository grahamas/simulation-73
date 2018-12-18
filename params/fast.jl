using Modeling, WCM, Meshes, Records, Simulating,
  CalculatedParameters, Analysis, WCMAnalysis,
  WCMConnectivity, WCMNonlinearity, WCMStimulus, WCMTarget
using WCM: WCMSpatial1D
using DifferentialEquations: Euler

NUM = Float64
if !(@isdefined UV)
  const UV = UnboundedVariable
  const BV = BoundedVariable
  const varying{T} = Union{T,BV{T}}
  const v = NUM
end
T= 3.0
simulation = Simulation(
        model = WCMSpatial1D(;
            pop_names = ["E", "I"],
            α = [1.1, 1.0],
            β = [1.1, 1.1],
            τ = [0.1, 0.18],
            space = Segment{v}(; n_points=201, extent=100),
            nonlinearity = pops(SigmoidNonlinearity{v}; a=[1.2, 1.0],
                                                        θ=[2.6, 8.0]),
            stimulus = add([
                            pops(SharpBumpStimulus{v}; strength=[1.2, 1.2],
                                                   window=[(0.5,0.65), (0.5,0.65)],
                                                   width=[2.81, 2.81]),
                            pops(GaussianNoiseStimulus{v,1}; SNR=[80.0, 80.0])]),
            connectivity = pops(ShollConnectivity{v};
                amplitude = [16.0 -18.2;
                             27.0 -4.0],
                spread = v[2.5 2.7;
                           2.7 2.5])
            ),
        solver = AutoSolver(;
            T = T
            #:dt => 0.02,
            ),
        analyses = Analyses{v}(;
          subsampler = SubSampler(;
               space_strides = [2],
               dt = 0.01
               ),
          plots = [
              # Animate(;
              #   fps = 20
              # )#,
              # NonlinearityPlot(;
              #   fn_bounds = (-1,15)
              # ),
              #SpaceTimePlot(),
              NeumanTravelingWavePlot(;
                dt = 0.1
              )
              ]
            ),
        output = SingleOutput(;
            root = "/home/grahams/Dropbox/simulation-73/results/",
            simulation_name = "fixing_neuman"
            )
        )

using JLD2

@save "parameters.jld2" simulation
