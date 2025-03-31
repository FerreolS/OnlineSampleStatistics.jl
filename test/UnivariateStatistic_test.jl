@testset "UnivariateStatistic" begin
    # Test for the UnivariateStatistic struct
    A = UnivariateStatistic(2, 0.5)
    B = UnivariateStatistic(1, Float32)

    @test A.rawmoments == [0.5, 0.0]
    @test A.weights == 1
    @test B.rawmoments == [0.0f0]
    @test B.weights == 0

    @test mean(A) == 0.5
    @test var(A) === NaN
    @test skewness(A) === NaN
    @test kurtosis(A) === NaN
    @test nobs(A) == 1

    C = zero(A)
    push!(C, ones(10))
    @test mean(C) == 1.0
    @test var(C) == 0.0

    D = zero(typeof(A))
    push!(D, zeros(10))

    merge!(C, D)
    @test mean(C) == 0.5
    @test var(C; corrected=false) == 0.25
end
