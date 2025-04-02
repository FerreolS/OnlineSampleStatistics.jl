@testset "UnivariateStatistic" begin
    # Test for the UnivariateStatistic struct
    A = UnivariateStatistic(2, 0.5)
    B = UnivariateStatistic(1, Float32)

    @test A == UnivariateStatistic(Float64, 2, 0.5)
    @test A.rawmoments == [0.5, 0.0]
    @test A.weights == 1
    @test B.rawmoments == [0.0f0]
    @test B.weights == 0

    @test @inferred mean(A) == 0.5
    @test @inferred isnan(var(A))
    @test @inferred isnan(skewness(A))
    @test @inferred isnan(kurtosis(A))
    @test @inferred nobs(A) == 1

    @testset "push! merge! tests for first two moments" begin
        C = zero(A)
        push!(C, ones(10))
        @test @inferred mean(C) == 1.0
        @test @inferred var(C) == 0.0

        D = zero(typeof(A))
        push!(D, zeros(10))

        @test UnivariateStatistic(Float64, 2, [1.0f0, 0.0f0]) == UnivariateStatistic(2, [1.0, 0.0])

        merge!(C, D)
        @test @inferred mean(C) == 0.5
        @test @inferred var(C; corrected=false) == 0.25

        E = UnivariateStatistic(2, Float32)
        push!(E, ones(Float32, 10))
        F = UnivariateStatistic(2, Float64)
        @inferred push!(F, zeros(Float32, 10))
        @inferred merge!(F, E)
        @test mean(F) == 0.5
        @test eltype(F) == Float64

        F = UnivariateStatistic(2, Float64)
        E = UnivariateStatistic(2, Float32)
        push!(F, zeros(Float32, 10))
        @test merge!(F, E) == UnivariateStatistic([0.0, 0.0], 10)

        F = UnivariateStatistic(1, Float64)
        E = UnivariateStatistic(1, Float32)
        push!(F, zeros(Float32, 10))
        push!(E, zeros(Float32, 10))
        @test @inferred merge(F, E) == UnivariateStatistic([0.0], 20)
        @test @inferred merge(E, F) == UnivariateStatistic([0.0], 20)
    end

    @testset "higher moments" begin
        x = 1e9 .+ (randn(10^6)) .^ 2
        A = UnivariateStatistic(4)
        push!(A, x)
        @test @inferred mean(A) ≈ mean(x)
        @test @inferred isapprox(var(A), var(x); rtol=1e-6)
        @test @inferred isapprox(var(A; corrected=false), var(x; corrected=false); rtol=1e-6)
        @test @inferred isapprox(skewness(A), skewness(x); rtol=1e-6)
        @test @inferred isapprox(kurtosis(A), kurtosis(x); rtol=1e-6)
    end
end
