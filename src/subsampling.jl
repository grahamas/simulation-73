
function scalar_to_idx_window(scalar_window::Tuple{T,T}, origin_idx::Int, max_idx::Int, Δscalar::T) where T
	unbounded_float_idx_window = origin_idx .+ scalar_window ./ Δscalar
	bounded_float_idx_window = max.(min.(unbounded_float_idx_window, max_idx), 1)
	return round.(Int, bounded_float_idx_window)
end

function scalar_to_idx_window(scalar_window::Tuple{T,T}, arr::AbstractArray{T,1}) where T
	(findfirst((scalar_window[1] .< arr) .| (scalar_window[1] .≈ arr)), findlast((scalar_window[2] .> arr) .| (scalar_window[2] .≈ arr)))
end

function subsampling_Δidx(Δsubsampled::T, Δsource::T) where T
	Δidx = max(1, round(Int, Δsubsampled / Δsource))
	return Δidx
end

function subsampling_idxs(Δsource::T, max_idx::Int; Δsubsampled::T=Δsource, scalar_window::Union{Tuple{T,T},Nothing}=nothing, origin_idx::Int=1) where T
	if scalar_window != nothing
		idx_window = scalar_to_idx_window(scalar_window, origin_idx, max_idx, Δsource)
	else
		idx_window = (origin_idx, max_idx)
	end
	idx_window[1]:subsampling_Δidx(Δsubsampled, Δsource):idx_window[2]
end

function subsampling_idxs(target::AbstractArray{T,1}, source::AbstractArray{T,1}) where T
	Δsource = source[2] - source[1]
	max_idx = length(source)

	Δsubsampled = target[2] - target[1]
	scalar_window = (target[1], target[end])
	idx_window = scalar_to_idx_window(scalar_window, source)

	Δidx = subsampling_Δidx(Δsubsampled, Δsource)

	idx_window[1]:Δidx:idx_window[2]
end

function subsampling_idxs(target::AbstractArray{T,1}, Δsource::T, max_idx::Int; origin_idx::Int=1) where T
	subsampling_idxs(Δsource, max_idx; scalar_window=(target[1], target[end]), Δsubsampled=target[2]-target[1], origin_idx=origin_idx)
end
function subsampling_time_idxs(t_target, solver)
    t_solver = time_span(solver)[1]:save_dt(solver):time_span(solver)[end]
    subsampling_idxs(t_target, t_solver)
end
function subsampling_space_idxs(x_target, model, solver)
    x_model = space_arr(model)[1:solver.space_save_every:end,1]
    subsampling_idxs(x_target, x_model)
end