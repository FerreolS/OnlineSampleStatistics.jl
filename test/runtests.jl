using OnlineSampleStatistic
using Test
using Aqua

@testset "OnlineSampleStatistic.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OnlineSampleStatistic)
    end
    # Write your tests here.
end
