using Test, Documenter, Simulation73
#DocMeta.setdocmeta!(Simulation73, :DocTestSetup, :(using Simulation73); recursive=true)
#doctest(Simulation73)

n_points_dx1 = 1001
extent_dx1 = 300.0
dx = extent_dx1 / n_points_dx1
mid_point = floor(Int, n_points_dx1 / 2) + (n_points_dx1 % 2)
circle_dx1 = PeriodicLattice{Float64,1}(; n_points=(n_points_dx1,), extent=(extent_dx1,))
line_dx1 = CompactLattice{Float64,1}(; n_points=(n_points_dx1,), extent=(extent_dx1,))
empty_circle_dx1 = zeros(size(circle_dx1)...)
empty_line_dx1 = zeros(size(line_dx1)...)

@testset "Stimulus" begin
    @testset "Non-stimulus" begin
        nostim = NoStimulusParameter{Float64}()
        nostim_action = nostim(circle_dx1)
        empty_circle_nostim_test = copy(empty_circle_dx1)
        nostim_action(empty_circle_nostim_test)
        @test all(empty_circle_nostim_test .== empty_circle_dx1)
    end
end
@testset "Populations" begin
    @testset "Initializing" begin
        single_pop = copy(empty_line_dx1)
        two_pops = population_repeat(single_pop, 2)
        derived_single_pop = population(two_pops, 1)
        @test all(size(derived_single_pop) .== size(single_pop))
    end
    @testset "Non-stimuli" begin
        nostim1 = NoStimulusParameter{Float64}()
        nostim2 = NoStimulusParameter{Float64}()
        nostim_pops = PopulationActionsParameters(nostim1, nostim2)
        nostim_actions = nostim_pops(line_dx1)
        single_pop = copy(empty_line_dx1)
        two_pops = population_repeat(single_pop, 2)
        @test nostim_actions(two_pops, two_pops, 1.0) != 0
    end
    @testset "Meta-parameters" begin
        nostim1 = NoStimulusParameter{Float64}()
        nostim2 = NoStimulusParameter{Float64}()
        nostim_pops = PopulationActionsParameters(nostim1, nostim2)
        nullstim_pops = NullifyParameter(nostim_pops)
        nullstim_actions = nullstim_pops(line_dx1)
        single_pop = copy(empty_line_dx1)
        two_pops = population_repeat(single_pop, 2)
        @test nullstim_actions(two_pops, two_pops, 0.0) != 0 
    end
end



