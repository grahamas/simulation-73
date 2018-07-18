module Targets

using CalculatedParameters

abstract type Target{T} <: Parameter{T} end
abstract type TargetFactory{T} <: Target{T} end
abstract type LossFunction{T} <: Target{T} end

export Target, TargetFactory, LossFunction

end