
#region AbstractModel Deconstruction & Reconstruction
#########################################################
######### MODEL DECONSTRUCTION & RECONSTRUCTION #########
# Deconstruct a model (which is naturally hierarchical)
# into a flat representation.
# Alter flat representation.
# Reconstruct hierarchical model.
#########################################################

"""
    set_deep_dx!(model_deconstruction, dxs, val)

Within an object deconstructed into an array, set one element to val,
where dxs indicate that element by successive indexing into nested arrays.
"""
function set_deep_dx!(model_deconstruction::Tuple{Type,Array}, dxs, val)
    tmp = model_deconstruction
    for dx in dxs[1:end-1]
        tmp = tmp[2][dx]
    end
    tmp[2][dxs[end]] = val
end

function deconstruct(v::AbstractVariable)
    val = default_value(v)
    return (typeof(val), val)
end

function deconstruct(val::V) where {V <: Union{Number, AbstractString}}
    return (typeof(val), val)
end

function deconstruct(m::M) where {M <: AbstractParameter}
    deconstruction = Tuple{Type,Any}[]
    for i_field in 1:nfields(m)
        substruct = deconstruct(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return (typeof(m), deconstruction)
end

function deconstruct(arr::AA) where {AA <: AbstractArray}
    deconstruction = map(deconstruct, arr)
    return (base_type(typeof(arr)), [deconstruction...]) # Remade in case static
end

function deconstruct(tup::Tuple)
    deconstruction = map(deconstruct, tup)
    return (base_type(typeof(tup)), [deconstruction...])
end

function var_deconstruct(val::AbstractVariable)
    return (typeof(val), val)
end

function var_deconstruct(val::Union{AbstractString,Number})
    return (typeof(val), val)
end

function var_deconstruct(m::M) where {M <: AbstractParameter}
    deconstruction = Tuple{Type,Any}[]
    for i_field in 1:nfields(m)
        substruct = var_deconstruct(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return (typeof(m), deconstruction)
end

function var_deconstruct(arr::AA) where {AA<:AbstractArray}
    deconstruction = map(var_deconstruct, arr)
    return (base_type(typeof(arr)), [deconstruction...]) # remade incase static
end

function var_deconstruct(tup::Tuple)
    deconstruction = map(var_deconstruct, tup)
    return (base_type(typeof(tup)), [deconstruction...])
end

function base_type(T::Int)
    T
end


function base_type(::Type{P}) where {P <: AbstractParameter}
    BTs = map(base_type, P.parameters)
    return (P.name.wrapper){BTs...}
end

function base_type(::Type{T}) where T <: Real
    return T
end

function base_type(::Type{V}) where {T, V<: Union{AbstractVariable{T}, T}}
    return T
end

function base_type(::Type{Array{T,N}}) where {T, N}
    BT = base_type(T)
    return Array{BT,N}
end

function base_type(::Type{Tuple{T,S}}) where {T,S}
    BT = base_type(T)
    BS = base_type(S)
    return Tuple{BT,BS}
end

function base_type(::Type{SA}) where {N,M,T,TUP,SA<:Union{SArray{TUP,T,N,M},MArray{TUP,T,N,M},SizedArray{TUP,T,N,M}}}
    BT = base_type(T)
    return SArray{TUP,BT,N,M}
end

function reconstruct(tup::Tuple{Type,<:Union{Number,AbstractString}})
    return tup[2]
end

function reconstruct(tup::Tuple{Type,Array})
    typ, arr = tup
    base_typ = base_type(typ)
    if typ <: Union{AbstractArray, Tuple}
        return base_typ(reconstruct.(arr))
    else
        return base_typ(reconstruct.(arr)...)
    end
end

############################################
##endregion
