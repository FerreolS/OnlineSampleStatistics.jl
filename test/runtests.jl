using OnlineSampleStatistics
using Test

@testset "OnlineSampleStatistics.jl" begin

    include("UnivariateStatistic_test.jl")
    include("IndependantStatistics_test.jl")
    include("show_test.jl")
end
