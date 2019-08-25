
Base.minimum(solution::DESolution) = minimum(map(minimum, solution.u))
Base.maximum(solution::DESolution) = maximum(map(maximum, solution.u))
