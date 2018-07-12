using Exploration, WC73, Meshes, Records, CalculatedParameters, WCMConnectivity, WCMNonlinearity, WCMStimulus, Targets
using WC73: WCMSpatial1D

if !isdefined(:UV)
  const UV = UnboundedVariable
  const BV = BoundedVariable
  const varying{T} = Union{T,BV{T}}
  const v = varying{Float64}
end

p_search = ParameterSearch(
        variable_model = WCMSpatial1D(;#{varying{Float64}}(;
            α = v[1.1, 1.0],
            β = v[1.1, 1.1],
            τ = v[0.1, 0.18],
            space = Segment{v}(; n_points=401, extent=100.5),
            nonlinearity = pops(SigmoidNonlinearity{v}; a=[1.2, 1.0], θ=[2.6, 8.0]),
            stimulus = pops(SharpBumpStimulus{v}; strength=[2.4,2.4], duration=[1.0,1.0], width=[3.0,3.0]),
            connectivity = pops(ShollConnectivity{v};
                amplitude = v[BV(16.0, (10.0,30.0)) -18.2;
                BV(27.0, (10.0,30.0)) -4.0],
                spread = v[2.5 2.7;
                2.7 2.5])
            ),
        solver = Solver(;
            T = 6,
            params = Dict(
                :dt => 0.001,
                #:alg_hints => [:stiff]
                )
            ),
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
        output = SingleOutput(;
            root = "/home/grahams/Dropbox/Research/simulation-73/results/",
            simulation_name = ""
            ),
        target_factory = DecayingWaveFactory(;
            decay = 2.5,
            speed = 10.0,
            timepoints=1:3,
            target_pop=1
            )
        )

using JLD

jldopen("parameters.jld", "w") do file
  addrequire.(file, [WC73, Meshes, Records, CalculatedParameters, WCMConnectivity, WCMNonlinearity, WCMStimulus])
  write(file, "p_search", p_search)
end
