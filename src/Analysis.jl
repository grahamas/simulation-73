module Analysis
module CalculatedParameters

# * Analysis types
struct Analyses{M <: Model}
	sampling::Sampling
	plots::Array{AbstractPlot}
end

struct Sampler{M <: Model}
	v_stride::Array{Integer}
end

abstract type AbstractPlot end

export analyse, Analyses

end
