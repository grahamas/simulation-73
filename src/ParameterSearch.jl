

function make_WilsonCowan73!(mesh, α, β, τ, nonlinearity_args, stimulus_args, connectivity_args) where T
    nonlinearity_fn = make_nonlinearity(mesh; nonlinearity_args...)
    stimulus_fn = make_stimulus(mesh; stimulus_args...)
    connectivity_fn = make_connectivity(mesh; connectivity_args...)
    function wc!(dA, A, p, t)
        @. dA = (-α * A + β * (1-A) * nonlinearity_fn(connectivity_fn(A) + stimulus_fn(t))) / τ
    end
end
