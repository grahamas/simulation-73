module Targets

using CalculatedParameters

abstract type TargetFactory{T} <: Parameter{T} end

export TargetFactory

end