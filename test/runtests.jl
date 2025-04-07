using OnlineSampleStatistics
using Test

@testset "OnlineSampleStatistics.jl" begin
    using OnlineStatsBase

    include("UnivariateStatistic_test.jl")
    include("IndependantStatistics_test.jl")
    # Write your tests here.
end
