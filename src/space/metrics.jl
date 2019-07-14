@doc """
    euclidean_metric(edge)

Return the distance between two points in euclidean space, given an edge between those points.

# Example
```jldoctest
julia> euclidean_metric( (5,1) )
4

julia> euclidean_metric( ((2,2), (5,-5)) )
(3, 7)
```
"""
euclidean_metric(edge::Tuple{T,T}) where T<:Number = abs(edge[1] - edge[2])
# FIXME should really be L2 norm
euclidean_metric(edge::Tuple{Tup,Tup}) where {T,N,Tup<:NTuple{N,T}} = abs.(edge[1] .- edge[2])



"""
    euclidean_metric_periodic(edge, period)

Return the distance between two points in euclidean space as in euclidean_metric, but let the space wrap with period.

# Example
```jldoctest
julia> euclidean_metric_periodic( (5,1), 3 )
3

julia> euclidean_metric_periodic( ((2,2), (5,-5)), (3,4) )
(0, 3)
```
"""
function euclidean_metric_periodic(edge::Tuple{T,T}, period::T) where T<:Number
    diff = euclidean_metric(edge)
    if diff > period / 2
        return period - diff
    else
        return diff
    end
end
function euclidean_metric_periodic(edge::Tuple{Tup,Tup}, periods::Tup) where {N,T,Tup<:NTuple{N,T}}
    diffs = euclidean_metric(edge)
    diffs = map(zip(diffs, periods)) do (diff, period)
        if diff > period / 2
            return period - diff
        else
            return diff
        end
    end
    return Tup(diffs)
end
