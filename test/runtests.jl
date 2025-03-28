using OnlineSampleStatistics
using Test
using Aqua

@testset "OnlineSampleStatistics.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OnlineSampleStatistics)
    end
    # Write your tests here.
end
