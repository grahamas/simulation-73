module Connectivity

# * Top factories

function make_connectivity_mx(mesh; name=error("Missing arg"), args...)
    connectivity_mx_factories = Dict(
    "sholl" => make_sholl_mx
    )
    return connectivity_mx_factories[name](mesh; args...)
end

# * Distance matrix
doc"""
This matrix contains values such that the $j^{th}$ column of the $i^{th}$ row
contains the distance between locations $i$ and $j$ in the 1D space dimension provided.
"""
function distance_matrix(xs::SpaceDim{DistT}) where {DistT <: Real}
    # aka Hankel, but that method isn't working in SpecialMatrices
    distance_mx = zeros(DistT, length(xs), length(xs))
    for i in range(1, length(xs))
        distance_mx[:, i] = abs.(xs - xs[i])
    end
    return distance_mx'
end
# * Sholl connectivity

function make_sholl_mx(mesh::FlatMesh; A, σ)
    space_dim = mesh.space.dims[1]
    dist_mx = distance_matrix(space_dim)
    step_size = step(space_dim)
    return sholl_matrix(A, σ, dist_mx, step_size)
end

@memoize function make_sholl_mx(, A, σ)

doc"""
We use an exponential connectivity function, inspired both by Sholl's
experimental work, and by certain theoretical considerations.

The interaction between two populations is entirely characterized by this
function and its two parameters: the amplitude (weight) and the spread
(σ). The spatial step size is also a factor, but as a computational concern
rather than a fundamental one.
"""
function sholl_matrix(amplitude::ValueT, spread::ValueT,
                      dist_mx::Interaction1DFlat{ValueT}, step_size::ValueT) where {ValueT <: Real}
    conn_mx = @. amplitude * step_size * exp(
        -abs(dist_mx / spread)
    ) / (2 * spread)
    return conn_mx
end
doc"""
This calculates a matrix of Sholl's exponential decay for each pair of
populations, thus describing all pairwise interactions. The result is a tensor
describing the effect of the source population at one location on the target
population in another location (indexed: `[tgt_loc, tgt_pop, src_loc,
src_pop]`). This works for arbitrarily many populations (untested) but only for
1D space.
"""
function sholl_connectivity(mesh::PopMesh{DistT}, W::InteractionParam{ValueT},
			    Σ::InteractionParam{ValueT})::InteractionTensor{ValueT} where {ValueT <: Real, DistT <: Real}
    xs = mesh.space.dims[1]
    N_x = length(xs)
    N_pop = size(W)[1]
    conn_tn = zeros(N_x, N_pop, N_x, N_pop)
    for tgt_pop in range(1,N_pop)
	for src_pop in range(1,N_pop)
	    conn_tn[:, tgt_pop, :, src_pop] .= sholl_matrix(W[tgt_pop, src_pop],
			  Σ[tgt_pop, src_pop], distance_matrix(xs), step(xs))
	end
    end
    return conn_tn
end
doc"""
In the two population case, flattening the tensor and using matrix
multiplication is 3x faster. This function provides exactly that.
"""
function sholl_connectivity(mesh::FlatMesh{ValueT}, args...) where {ValueT <: Real}
    # Why didn't I provide an unflattened mesh in the first place?
    sholl_connectivity(unflatten(mesh), args...) |> flatten_sholl
end
function flatten_sholl(tensor)::Interaction1DFlat
    N_x, N_p = size(tensor)[1:2]
    @assert N_p < N_x
    @assert size(tensor) == (N_x, N_p, N_x, N_p)
    flat = zeros(eltype(tensor), N_x*N_p, N_x*N_p)
    for i in 1:N_p
        for j in 1:N_p
            flat[(1:N_x)+((i-1)*N_x), (1:N_x)+((j-1)*N_x)] = tensor[:,i,:,j]
        end
    end
    return flat
end

end
