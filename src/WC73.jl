module WC73

dbg(x) = @schedule println(x)
dbg(x::Array) = @schedule println("type $(typeof(x)), size $(size(x))")

include("Calculated.jl")
include("SimulationTypes.jl")

include("Connectivity.jl")
include("Stimulus.jl")
include("Nonlinearity.jl")

include("Parameters.jl")
include("Simulation.jl")
include("Analysis.jl")

end
