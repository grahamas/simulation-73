module WC73

dbg(x) = @schedule println(x)
dbg(x::Array) = @schedule println("type $(typeof(x)), size $(size(x))")

include("SimulationTypes.jl")
include("Parameters.jl")
include("Simulation.jl")
include("Analysis.jl")

end
