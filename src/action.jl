
abstract type AbstractAction{T} end
abstract type AbstractInteraction{T} <: AbstractAction{T} end
abstract type AbstractSpaceAction{T,N} end
abstract type AbstractSpaceInteraction{T,N} <: AbstractSpaceAction{T,N} end
abstract type AbstractPopulationP{NP,AP} end

struct NullAction{T} <: AbstractAction{T} end
(na::NullAction)(args...) = nothing

struct NullifyParameter{T} <: AbstractParameter{T} end
NullifyParameter(param::AbstractParameter{T}) where T = NullifyParameter{T}()
(np::NullifyParameter{T})(space::AbstractSpace{T}) where {T} = NullAction{T}()
