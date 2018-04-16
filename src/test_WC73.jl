using Base.Test
push!(LOAD_PATH, ".")
using WC73
import WC73: simple_sigmoid_fn, sigmoid_fn
@testset "Sigmoids" begin
    @test simple_sigmoid_fn(0,1,0) == 0.5
    @test simple_sigmoid_fn(0,1,0.5) ≈ 0.37754066879814
    @test simple_sigmoid_fn(1,1,1) == 0.5
    @test sigmoid_fn(0,1,0) == 0.0
    @test sigmoid_fn(0,1,0.5) == 0.0
    @test sigmoid_fn(1,1,1) ≈ 0.231058578630049
end
@testset "Stimulus" begin
    @test_skip true
end
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
