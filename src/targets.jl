module Targets

using CalculatedParameters
using Modeling
using Meshes

abstract type Target{T} <: Parameter{T} end
abstract type TargetFactory{T} <: Target{T} end
abstract type LossFunction{T} <: Target{T} end

function loss(factory::TargetFactory, model::Model)
    target_data = factory(model)
    weights = zeros(target_data)
    weights[Calculated(model.space).value .>= 0,factory.target_pop,:] .= 1
    L2Loss(factory.timepoints, target_data; data_weight=weights)
end

export Target, TargetFactory, LossFunction, loss

end