using Modeling, WCM, Meshes, Records, Simulating,
  CalculatedParameters, Analysis, WCMAnalysis,
  WCMConnectivity, WCMNonlinearity, WCMStimulus, WCMTarget
using WCM: WCMSpatial1D

if !(@isdefined UV)
  const UV = UnboundedVariable
  const BV = BoundedVariable
  const varying{T} = Union{T,BV{T}}
  const v = Float64
end
T= 0.5
M = WCMSpatial1D
simulation = Simulation{M}(
        model = M(;
            pop_names = ["E", "I"],
            α = [1.1, 1.0],
            β = [1.1, 1.1],
            τ = [0.1, 0.18],
            P = [0.0, 0.0],
            space = Segment{v}(; n_points=281, extent=100),
            nonlinearity = pops(SigmoidNonlinearity{v}; a=[1.2, 1.0],
                                                        θ=[2.6, 8.0]),
            stimulus = pops(SharpBumpStimulus{v}; strength=[1.2, 1.2],
                                                  duration=[0.55, 0.55],
                                                  width=[2.81, 2.81]),
            connectivity = pops(ShollConnectivity{v};
                amplitude = [16.0 -18.2;
                             27.0 -4.0],
                spread = v[2.5 2.7;
                           2.7 2.5])
            ),
        solver = Solver(;
            T = T,
            params = Dict(
                :dt => 0.1,
                #:dense => true
                #:alg_hints => [:stiff]
                )
            ),
        analyses = Analyses{M}(;
          subsampler = SubSampler{M}(;
               space_strides = [4],
               dt = 0.05
               ),
          plots = [
              Animate(;
                fps = 20
                ),
              NonlinearityPlot(;
                fn_bounds = (-1,15)),
              SpaceTimePlot()
              ]
            ),
        output = SingleOutput(;
            root = "/home/grahams/Dropbox/simulation-73/results/",
            simulation_name = "ztmp_fast"
            )
        )

using JLD2

@save "parameters.jld2" simulation
