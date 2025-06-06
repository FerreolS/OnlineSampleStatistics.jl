@testset "UnivariateStatistic" begin
    # Test for the UnivariateStatistic struct
    @testset "UnivariateStatistic Constructors" begin

        # Test constructor with default type
        A = UnivariateStatistic(0.5, 2)
        @test typeof(A) <: UnivariateStatistic
        @test @inferred isnan(var(A))
        @test @inferred nobs(A) == 1

        # Test constructor with custom type
        B = UnivariateStatistic(Float32, 1)
        @test typeof(B) <: UnivariateStatistic


        @test A == UnivariateStatistic(Float64, 0.5, 2)
        @test A.rawmoments == [0.5, 0.0]
        @test A.weights == 1
        @test B.rawmoments == [0.0f0]
        @test B.weights == 0

        A = UnivariateStatistic([0.0, 1, 2], 2)
        @test nobs(A) == 3
        @test @inferred mean(A) == 1.0
    end


    @testset "OnlineStatsBase API" begin
        using OnlineStatsBase

        A = UnivariateStatistic(1, 2)
        @test value(A) == [1, 0.0]
        @test empty!(A) == zero(A)
        fit!(A, ones(10))
        @test A == UnivariateStatistic(10, [1.0, 0.0])
    end

    @testset "fit! merge! tests for first two moments" begin
        A = UnivariateStatistic(0.5, 2)
        C = zero(A)
        fit!(C, ones(10))
        @test @inferred mean(C) == 1.0
        @test @inferred var(C) == 0.0
        @test @inferred std(C) == 0.0
        @test_throws ArgumentError UnivariateStatistic(-1, [1.0])


        D = zero(typeof(A))
        fit!(D, zeros(10))

        @test UnivariateStatistic(Float64, [1.0f0, 0.0f0], 2) == UnivariateStatistic([1.0, 0.0], 2)

        merge!(C, D)
        @test @inferred mean(C) == 0.5
        @test @inferred var(C; corrected=false) == 0.25

        E = UnivariateStatistic(Float32, 2)
        fit!(E, ones(Float32, 10))
        F = UnivariateStatistic(Float64, 2)
        @inferred fit!(F, zeros(Float32, 10))
        @inferred merge!(F, E)
        @test mean(F) == 0.5
        @test eltype(F) == Float64

        F = UnivariateStatistic(Float64, 2)
        E = UnivariateStatistic(Float32, 2)
        fit!(F, zeros(Float32, 10))
        @test merge!(F, E) == UnivariateStatistic(10, [0.0, 0.0])

        F = UnivariateStatistic(Float64, 1)
        E = UnivariateStatistic(Float32, 1)
        fit!(F, zeros(Float32, 10))
        fit!(E, zeros(Float32, 10))
        @test @inferred merge(F, E) == UnivariateStatistic(20, [0.0])
        @test @inferred merge(E, F) == UnivariateStatistic(20, [0.0])

    end

    @testset "higher moments" begin
        x = 1e9 .+ (randn(10^6)) .^ 2
        A = UnivariateStatistic(4)
        @test @inferred isnan(skewness(A))
        @test @inferred isnan(kurtosis(A))
        fit!(A, x)
        @test @inferred(order(A)) == 4
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

    @testset "Complex numbers" begin
        x = Complex.(1e9 .+ (randn(10^6)) .^ 2, randn(10^6))
        @test_throws ArgumentError UnivariateStatistic(x, 4)

        A = UnivariateStatistic(x, 2)
        x1 = x[1:10^5]
        x2 = x[(10^5+1):end]
        B = UnivariateStatistic(x1, 2)
        C = UnivariateStatistic(x2, 2)
        merge!(B, C)
        @test isapprox(A.rawmoments, B.rawmoments; rtol=1e-6)
    end

    @testset "error handling" begin
        @test_throws ArgumentError UnivariateStatistic(2, -1)
        @test_throws ArgumentError UnivariateStatistic(2, 0)
        @test_throws ArgumentError UnivariateStatistic(2, 1, -1)
        @test_throws ArgumentError UnivariateStatistic(2, 1, 0)
    end

    @testset "Transducers Tests" begin
        using Transducers
        x = 1e9 .+ (randn(10^6)) .^ 2
        A = UnivariateStatistic(4)
        fit!(A, x)
        B = foldxt(UnivariateStatistic(4), x)
        @test isapprox(A, B; rtol=1e-6)

        A = fit!(UnivariateStatistic(Complex{Float64}, 2), x)
        B = foldxt(UnivariateStatistic(Complex{Float64}, 2), Complex.(x))
        @test isapprox(A, B; rtol=1e-6)

    end

    @testset "Weighted Data" begin
        UnivariateStatistic(ones(Int, 10), (randn(10) .> 0), 4)
        x = randn(1_000)
        w = rand(size(x)...)
        A = UnivariateStatistic(x, w, 4)
        @test nobs(A) ≈ sum(w)
        @test mean(A) ≈ sum(w .* x) ./ sum(w)
        @test isapprox(var(A, corrected=false), sum(w .* (x .- mean(A)) .^ 2) ./ sum(w); rtol=1e-6)
        @test_logs (:warn, "The number of samples is not an integer. The variance is not corrected.") var(A)


        A = UnivariateStatistic(Float64, Float64, 1)
        fit!(A, x, 2.0)
        @test @inferred(nobs(A)) ≈ 2 * length(x)
        @test @inferred(mean(A)) ≈ mean(x)

        A = UnivariateStatistic(Float64, Float64, 2)
        fit!(A, Float32.(x), 2.0)
        @test @inferred(nobs(A)) ≈ 2 * length(x)
        @test @inferred(mean(A)) ≈ mean(Float32.(x))
        @test @inferred(var(A, corrected=false)) ≈ var(x, corrected=false)

    end
    @testset "Weighted Data with Transducers" begin
        using Transducers
        x = randn(1_000)
        w = rand(size(x)...)
        A = @inferred(fit!(UnivariateStatistic(Float64, Float64, 4), zip(x, w)))
        B = foldxt(UnivariateStatistic(Float64, Float64, 4), zip(x, w))
        @test A ≈ B
    end

end
