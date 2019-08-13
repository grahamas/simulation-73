
Base.minimum(solution::DESolution) = minimum(map(minimum, solution.u))
Base.maximum(solution::DESolution) = maximum(map(maximum, solution.u))

"""
    pop_frame(solution, pop_dx, time_dx)

Return spatial frame for a given population `pop_dx` and time `time_dx`.
"""
@generated function pop_frame(solution::ODESolution{T,NPT,<:AbstractArray{<:AbstractArray{T,NP},1}}, pop_dx::Int, time_dx::Int) where {T,NP,NPT}
    N = NP - 1
    colons = [:(:) for i in 1:N]
    :(solution[pop_dx, $(colons...), time_dx])
end
