
Base.minimum(solution::DESolution) = minimum(map(minimum, solution.u))
Base.maximum(solution::DESolution) = maximum(map(maximum, solution.u))

"""
    pop_frame(solution, pop_dx, time_dx)

Return spatial frame for a given population `pop_dx` and time `time_dx`.
"""
@generated function pop_frame(solution::DiffEqBase.AbstractODESolution{T,NPT}, pop_dx::Int, time_dx::Int) where {T,NPT}
    N = NPT - 2 # N + pop + time
    colons = [:(:) for i in 1:N]
    :(solution[pop_dx, $(colons...), time_dx])
end
