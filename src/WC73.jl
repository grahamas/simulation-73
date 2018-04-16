module WC73
# * Load Modules
using TensorOperations

using DifferentialEquations

using Base.Test

import Base.Dates

# Add repo path to load path and load local packages.
push!(LOAD_PATH, joinpath(ENV["HOME"], "gits", "simulation-73"))
using SimulationTypes, WC73Parameters

# * Simulation function definitions
# ** Sigmoid functions
doc"""
The sigmoid function is defined
```math
\begin{align}
\mathcal{S}(x) = \frac{1}{1 + \exp(-a(x - \theta))}
\end{align}
```
where $a$ describes the slope's steepness and $\theta$ describes translation of the slope's center away from zero.

This is "simple" because in practice we use the rectified sigmoid.
"""
function simple_sigmoid_fn(x, a, theta)
    return @. (1 / (1 + exp(-a * (x - theta))))
end

doc"""
A rectified version of `simple_sigmoid_fn`.

In practice, we use rectified sigmoid functions because firing rates cannot be negative.

TODO: Rename to rectified_sigmoid_fn.
"""
function sigmoid_fn(x, a, theta)
    return max.(0, simple_sigmoid_fn(x, a, theta) .- simple_sigmoid_fn(0, a, theta))
end

# ** Stimulus functions
doc"""
    make_stimulus_fn(mesh; name, stimulus_args...)

A factory taking the domain (`mesh`) and `name` of a stimulus and returning the function
defined to be associated with that name mapped over the given domain.
"""
function make_stimulus_fn(mesh; name=nothing, stimulus_args...)
    stimulus_factories = Dict(
        "smooth_bump" => smooth_bump_factory,
        "sharp_bump" => sharp_bump_factory
    )
    return stimulus_factories[name](mesh; args...)
end

# *** Smooth bump
"Implementation of smooth_bump_frame used in smooth_bump_factory."
function make_smooth_bump_frame(mesh_coords::Array{DistType}, width::DistType, strength::NumType, steepness::NumType)
    @. strength * (simple_sigmoid_fn(mesh_coords, steepness, -width/2) - simple_sigmoid_fn(mesh_coords, steepness, width/2))
end

"""
The smooth bump is a smooth approximation of the sharp impulse defined
elsewhere. It is smooth in both time and space. It is constructed essentially
from three sigmoids: Two coplanar in space, and one orthogonal to those in
time. The two in space describe a bump: up one sigmoid, then down a negative
sigmoid. The one in time describes the decay of that bump.

This stimulus has the advantages of being 1) differentiable, and 2) more
realistic. The differentiabiilty may be useful for the automatic solvers that
Julia has, which can try to automatically differentiate the mutation function
in order to improve the solving.
"""
function smooth_bump_factory(mesh::AbstractMesh;
                             width=nothing, strength=nothing, duration=nothing,
                             steepness=nothing)
    # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_smooth_bump_frame(coords(mesh), width, strength, steepness)
    return (t) -> @. on_frame * (1 - simple_sigmoid_fn(t, steepness, duration))
end

# *** Sharp bump
# TODO Understand these functions again.....
"Implementation of sharp_bump_frame used in sharp_bump_factory"
function make_sharp_bump_frame(mesh::PopMesh, width::DistType, strength::NumType)
    mesh_coords = coords(mesh)
    frame = zeros(mesh_coords)
    mid_point = 0     # half length, half width
    half_width = width / 2      # using truncated division
    xs = mesh_coords[:,1]   # Assumes all pops have same mesh_coords
    start_dx = find(xs .>= mid_point - half_width)[1]
    stop_dx = find(xs .<= mid_point + half_width)[end]
    frame[start_dx:stop_dx,:] = strength
    return frame
end
function make_sharp_bump_frame(mesh::FlatMesh, args...)
    structured_frame = make_sharp_bump_frame(mesh.pop_mesh, args...)
    flat_frame = structured_frame[:] # Works because FlatMesh must have 1D PopMesh
    return flat_frame
end
"""
The "sharp bump" is the usual theoretical impulse: Binary in both time and
space. On, then off.
"""
function sharp_bump_factory(mesh; width=nothing, strength=nothing, duration=nothing)
        # WARNING: Defaults are ugly; Remove when possible.
    on_frame = make_sharp_bump_frame(mesh, width, strength)
    off_frame = zeros(on_frame)
    return (t) -> (t <= duration) ? on_frame : off_frame
end
# *** Stimulus testing
@testset "Stimulus" begin
    @test_skip true
end

# ** Connectivity functions

doc"""
This matrix contains values such that the $j^{th}$ column of the $i^{th}$ row
contains the distance between locations $i$ and $j$ in the 1D space dimension provided.
"""
function distance_matrix(xs::SpaceDim)
    # aka Hankel, but that method isn't working in SpecialMatrices
    distance_mx = zeros(eltype(xs), length(xs), length(xs))
    for i in range(1, length(xs))
        distance_mx[:, i] = abs.(xs - xs[i])
    end
    return distance_mx'
end
doc"""
We use an exponential connectivity function, inspired both by Sholl's
experimental work, and by certain theoretical considerations.

The interaction between two populations is entirely characterized by this
function and its two parameters: the amplitude (weight) and the spread
(σ). The spatial step size is also a factor, but as a computational concern
rather than a fundamental one.
"""
function sholl_matrix(amplitude::NumType, spread::NumType,
                      dist_mx::Array{NumType,2}, step_size::NumType)
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
function sholl_connectivity(mesh::PopMesh, W::Array{NumType,2},
			    Σ::Array{NumType,2})::InteractionTensor
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
function sholl_connectivity(mesh::FlatMesh, args...)
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

# *** Connectivity testing
@testset "Connectivity" begin
    @testset "Distance Matrix" begin
        @test_skip true
    end
    import WC73: sholl_matrix, distance_matrix
    @testset "Sholl Matrix" begin
        xs = linspace(-1.0,1.0,3)
        @test all(.≈(sholl_matrix(1.0, 1.0, distance_matrix(xs), step(xs)), [0.5         0.18393972  0.06766764;
                                                   0.18393972  0.5         0.18393972;
                                                   0.06766764  0.18393972  0.5       ], atol=1e-6))
    end
     import WC73: sholl_connectivity, PopMesh, flatten
     @testset "Sholl tensor" begin
           weights = [1.0 2.0; 3.0 4.0]
           spreads = [0.1 0.2; 0.3 0.4]
           mesh = PopMesh([Dict(:N => 3, :extent => 2)], 2)
           observed = sholl_connectivity(flatten(mesh), weights, spreads)
           expected =      [  5.00000000e+00   2.26999649e-04   1.03057681e-08   5.00000000e+00   3.36897350e-02   2.26999649e-04 ;
     2.26999649e-04   5.00000000e+00   2.26999649e-04   3.36897350e-02   5.00000000e+00   3.36897350e-02 ;
     1.03057681e-08   2.26999649e-04   5.00000000e+00   2.26999649e-04   3.36897350e-02   5.00000000e+00 ;
    5.          0.17836997  0.00636317  5.          0.41042499  0.03368973 ;
    0.17836997  5.          0.17836997  0.41042499  5.          0.41042499 ;
    0.00636317  0.17836997  5.          0.03368973  0.41042499  5.         ]
           println(observed)
           @test all(.≈(observed, expected, atol=1e-6))
     end
end

# ** Core differential functions
#=
Two versions of the differential equation function, both implementing the same
equation described in the Equations section. One is general, putatively
operating on any number of populations (though must still be 1D!). It is much
less efficient in the 2-pop case than the latter, which takes a matrix instead
of a tensor as the interaction term.
=#
WilsonCowan73!(dA,A::SpaceState1D,p::WilsonCowan73Params{InteractionTensor},t) = begin
    # Use dA as intermediate variable for tensor op since it is preallocated
    @tensor dA[x_tgt, pop_tgt] = p.W[x_tgt, pop_tgt, x_src, pop_src] * A[x_src, pop_src]
    dA .= (-p.α .* A + p.β .* (1 .- A) .* sigmoid_fn(dA + p.stimulus_fn(t), p.a, p.θ)) ./ p.τ
end

doc"""
This simulates a single timestep of the 1D spatial equation
```math
\begin{align}
\tau_E \partial_t E(x,t) &= -\alpha_E E(x,t) + \beta_E (1 - E(x,t)) \cS_E \left( W_{EE}(X) \conv E(x,t) + W_{EI}(X) \conv I(x,t) + P_E(x,t)\right)\\
\tau_I \partial_t I(x,t) &= -\alpha_I I(x,t) + \beta_I (1 - I(x,t)) \cS_I \left( W_{IE}(X) \conv E(x,t) +  W_{II}(X) \conv I(x,t) + P_I(x,t)\right)
\end{align}
```
"""
WilsonCowan73!(dA,A::SpaceState1DFlat,p::WilsonCowan73Params{Interaction1DFlat},t) = begin
    # Use dA as intermediate variable for tensor op since it is preallocated
    dA .= (-p.α .* A + p.β .* (1 .- A) .* sigmoid_fn(p.W*A + p.stimulus_fn(t), p.a, p.θ)) ./ p.τ
end
# * Solver function definition

function solve_WilsonCowan73(; model=nothing, solver=nothing, other...)
    solver_params = deepcopy(Dict{Any,Any}(solver))
    model_params = WilsonCowan73Params(model)

    u0 = zeros(model_params.mesh)
    tspan = (0.0, pop!(solver_params, :T))
    prob::DEProblem = ODEProblem(WilsonCowan73!, u0, tspan, model_params)

    if :dt in keys(solver_params)
        solver_params[:alg] = Euler()
        solver_params[:adaptive] = false
    elseif :stiff in keys(solver_params)
        if pop!(solver_params, :stiff) > 0
            solver_params[:alg_hints] = [:stiff]
        end
    end

    if :stiff in keys(solver_params)
        error("Incompatible solver parameters.")
    end

    soln::DESolution = solve(prob; solver_params...)

    return soln
end

# * Export and END
export solve_WilsonCowan73
end
