@testset "UnivariateStatistic" begin
    # Test for the UnivariateStatistic struct
    A = UnivariateStatistic(0.5, 2,)
    B = UnivariateStatistic(Float32, 1)

    @test A == UnivariateStatistic(Float64, 0.5, 2)
    @test A.rawmoments == [0.5, 0.0]
    @test A.weights == 1
    @test B.rawmoments == [0.0f0]
    @test B.weights == 0

    @test @inferred mean(A) == 0.5
    @test @inferred isnan(var(A))
    @test @inferred nobs(A) == 1

    @testset "push! merge! tests for first two moments" begin
        C = zero(A)
        push!(C, ones(10))
        @test @inferred mean(C) == 1.0
        @test @inferred var(C) == 0.0

        D = zero(typeof(A))
        push!(D, zeros(10))

        @test UnivariateStatistic(Float64, [1.0f0, 0.0f0], 2) == UnivariateStatistic([1.0, 0.0], 2)

        merge!(C, D)
        @test @inferred mean(C) == 0.5
        @test @inferred var(C; corrected=false) == 0.25

        E = UnivariateStatistic(Float32, 2)
        push!(E, ones(Float32, 10))
        F = UnivariateStatistic(Float64, 2)
        @inferred push!(F, zeros(Float32, 10))
        @inferred merge!(F, E)
        @test mean(F) == 0.5
        @test eltype(F) == Float64

        F = UnivariateStatistic(Float64, 2)
        E = UnivariateStatistic(Float32, 2)
        push!(F, zeros(Float32, 10))
        @test merge!(F, E) == UnivariateStatistic(10, [0.0, 0.0])

        F = UnivariateStatistic(Float64, 1)
        E = UnivariateStatistic(Float32, 1)
        push!(F, zeros(Float32, 10))
        push!(E, zeros(Float32, 10))
        @test @inferred merge(F, E) == UnivariateStatistic(20, [0.0])
        @test @inferred merge(E, F) == UnivariateStatistic(20, [0.0])
    end

    @testset "higher moments" begin
        x = 1e9 .+ (randn(10^6)) .^ 2
        A = UnivariateStatistic(4)
        @test @inferred isnan(skewness(A))
        @test @inferred isnan(kurtosis(A))
        push!(A, x)
        @test @inferred mean(A) ≈ mean(x)
        @test @inferred isapprox(var(A), var(x); rtol=1e-6)
        @test @inferred isapprox(var(A; corrected=false), var(x; corrected=false); rtol=1e-6)
        @test @inferred isapprox(skewness(A), skewness(x); rtol=1e-6)
        @test @inferred isapprox(kurtosis(A), kurtosis(x); rtol=1e-6)

        x1 = x[1:10^5]
        x2 = x[(10^5+1):end]
        B = UnivariateStatistic(x1, 4)
        C = UnivariateStatistic(x2, 4)
        merge!(B, C)
        @test isapprox(A.rawmoments, B.rawmoments; rtol=1e-6)
    end
end
