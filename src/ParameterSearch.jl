
module ParameterSearch

struct CalculatedWilsonCowan73{C<:Connectivity,
                               N<:Nonlinearity,S<:Stimulus}
    connectivity::Calculated{C}
    nonlinearity::Calculated{N}
    stimulus::Calculated{S}
end
function CalculatedWilsonCowan73(space::Space, connectivity::C,
                                 nonlinearity::N, stimulus::S) where {C<:Connectivity,
                                                                      N<:Nonlinearity, S<:Stimulus}
    CalculatedWilsonCowan73{C,N,T}(
        calculate(connectivity, space),
        calculate(nonlinearity),
        calculate(stimulus,space))
end

function update!(cwc::CalculatedWilsonCowan73{C,N,S}, p)
    connectivity = C(p)
    nonlinearity = N(p)
    stimulus = N(p)
    return any([update!(cwc, connectivity),
                update!(cwc, nonlinearity),
                update!(cwc, stimulus)])
end
update!(cwc::CalculatedWilsonCowan73, c::Connectivity) = update!(cwc.connectivity, c)
update!(cwc::CalculatedWilsonCowan73, n::Nonlinearity) = update!(cwc.nonlinearity, n)
update!(cwc::CalculatedWilsonCowan73, s::Stimulus) = update!(cwc.stimulus, s)
# No space changes allowed

function make_WilsonCowan73!(mesh, α, β, τ, nonlinearity_args, stimulus_args, connectivity_args) where T
    nonlinearity_fn = make_nonlinearity(mesh; nonlinearity_args...)
    stimulus_fn = make_stimulus(mesh; stimulus_args...)
    connectivity_fn = make_connectivity(mesh; connectivity_args...)
    function wc!(dA, A, p, t)
        @. dA = (-α * A + β * (1-A) * nonlinearity_fn(connectivity_fn(A) + stimulus_fn(t))) / τ
    end
end

function make_make_WilsonCowan73(space::Space, connectivity::Connectivity,
                                 nonlinearity::Nonlinearity, stimulus::Stimulus)
    cwc = CalculatedWilsonCowan73(space, connectivity, nonlinearity, stimulus)
    function make_WilsonCowan73(prob::ODEProblem, p)
        τ = p.τ
        α = p.α
        β = p.β
        update!(calculated_connectivity, connectivity)
        update!(calculated_nonlinearity, nonlinearity)
        update!(calculated_stimulus, stimulus)
        function WilsonCowan73!(dA, A, p, t)
            @. dA = (-α * A + β * (1-A) * nonlinearity_fn(connectivity_fn(A) + stimulus_fn(t))) / τ
        end
    end
end
end
