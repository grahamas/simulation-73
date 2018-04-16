@testset "Sigmoids" begin
    using WC73: simple_sigmoid_fn, sigmoid_fn
    @test simple_sigmoid_fn(0,1,0) == 0.5
    @test simple_sigmoid_fn(0,1,0.5) ≈ 0.37754066879814
    @test simple_sigmoid_fn(1,1,1) == 0.5
    @test sigmoid_fn(0,1,0) == 0.0
    @test sigmoid_fn(0,1,0.5) == 0.0
    @test sigmoid_fn(1,1,1) ≈ 0.231058578630049
end
