using ZippedArrays, StructuredArrays

@testset "IndependentStatistic Tests" begin
    # Test constructor with default type
    @testset "Constructor Tests" begin
        A = IndependentStatistic(5, (3, 3))
        @test size(A) == (3, 3)
        @test typeof(A) <: ZippedArray
    end

    # Test constructor with custom type
    @testset "other Type Constructor" begin
        A = IndependentStatistic(Float32, 2, (4, 4))
        @test size(A) == (4, 4)
        @test typeof(A) <: ZippedArray
    end

    @testset "Mean & variance Calculation" begin
        using OnlineSampleStatistics: get_rawmoments, weights, get_moments
        data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        A = IndependentStatistic(1, data)
        @test @inferred(mean(A)) ≈ mean.(A)

        A = IndependentStatistic(5, (3, 3))
        fit!(A, data)
        @test A == IndependentStatistic(5, data)
        @test @inferred(order(A)) == 5
        B = [UnivariateStatistic(5, d) for d in data]
        @test @inferred(order(B)) == 5
        @test A == B
        @test @inferred(get_rawmoments(A, 1)) == @inferred(get_rawmoments(B, 1))
        @test @inferred(get_rawmoments(A)) == @inferred(get_rawmoments(B))
        @test @inferred(get_moments(A, 1)) == @inferred(get_moments(B, 1))
        @test @inferred(weights(A)) == @inferred(weights(B))
        @test @inferred(nobs(A)) == @inferred(nobs(B))
        @test @inferred(mean(A)) == data

        @inferred(fit!(A, -1 .* data))
        @test @inferred(var(A; corrected = false)) == data .^ 2
        @test @inferred(var(A)) == 2 .* data .^ 2
        @test @inferred(skewness(A)) == zeros(Float64, 3, 3)
        @test @inferred(kurtosis(A)) == -2.0 * ones(3, 3)
        @test @inferred(var(A)) ≈ var.(A)
        @test @inferred(kurtosis(A)) ≈ kurtosis.(A)
        @test @inferred(skewness(A)) ≈ skewness.(A)
        @test @inferred(mean(A)) ≈ mean.(A)
        @test @inferred(var(A; corrected = false)) ≈ var.(A; corrected = false)

        C = IndependentStatistic(5, data, trues(size(data)...))
        @test @inferred(mean(C)) ≈ data
    end

    @testset "weightless moments estimation" begin
        sz = (2, 3)
        x = randn(sz..., 10)
        A = IndependentStatistic(1, sz)
        B = IndependentStatistic(Float64, 1, Int, sz)
        fit!(A, x)
        fit!(B, x)
        @test mean(A) ≈ dropdims(mean(x; dims = 3), dims = 3)
        @test mean(B) ≈ dropdims(mean(x; dims = 3), dims = 3)

        A = IndependentStatistic(2, sz)
        B = IndependentStatistic(Float64, 2, Int, sz)
        fit!(A, x)
        fit!(B, x)
        @test mean(A) ≈ dropdims(mean(x; dims = 3), dims = 3)
        @test mean(B) ≈ dropdims(mean(x; dims = 3), dims = 3)
        @test var(A) ≈ dropdims(var(x; dims = 3), dims = 3)
        @test var(B) ≈ dropdims(var(x; dims = 3), dims = 3)

        A = IndependentStatistic(4, sz)
        B = IndependentStatistic(Float64, 4, Int, sz)
        fit!(A, x)
        fit!(B, x)
        @test mean(A) ≈ dropdims(mean(x; dims = 3), dims = 3)
        @test mean(B) ≈ dropdims(mean(x; dims = 3), dims = 3)
        @test var(A) ≈ dropdims(var(x; dims = 3), dims = 3)
        @test var(B) ≈ dropdims(var(x; dims = 3), dims = 3)
        @test skewness(A) ≈ dropdims(mapslices(skewness, x; dims = 3), dims = 3)
        @test skewness(B) ≈ dropdims(mapslices(skewness, x; dims = 3), dims = 3)
        @test kurtosis(A) ≈ dropdims(mapslices(kurtosis, x; dims = 3), dims = 3)
        @test kurtosis(B) ≈ dropdims(mapslices(kurtosis, x; dims = 3), dims = 3)


    end
    # Test error handling for invalid dimensions
    @testset "Error Handling" begin
        A = IndependentStatistic(2, (3, 3))
        @test_throws DimensionMismatch fit!(A, [1.0 2.0; 3.0 4.0])
    end

    @testset "Weighted Data" begin
        x = randn(2, 3, 10)
        A = IndependentStatistic(2, x, trues(size(x)...))
        B = IndependentStatistic(2, x)
        @test @inferred(weights(A)) == @inferred(weights(B))
        @test @inferred(nobs(A)) == @inferred(nobs(B))
        @test @inferred(mean(A)) == mean(B)


        dims = 3
        A = IndependentStatistic(2, x, trues(size(x)...); dims = dims)
        B = IndependentStatistic(2, x; dims = dims)
        @test @inferred(weights(A)) == @inferred(weights(B))
        @test @inferred(nobs(A)) == @inferred(nobs(B))
        @test @inferred(mean(A)) == mean(B)
        @test @inferred(var(A, corrected = false)) ≈ var(B, corrected = false)


        w = (rand(size(x)...) .> 0.1)
        A = IndependentStatistic(2, x, w; dims = dims)

        @test mean(A) ≈ sum(w .* x; dims = dims) ./ sum(w; dims = dims)
        @test var(A; corrected = false) ≈ sum(w .* (x .- mean(A)) .^ 2; dims = dims) ./ sum(w; dims = dims)

        B = IndependentStatistic(Float64, 2, Int, size(A)[1:2])
        fit!(B, x, w)
        @test @inferred(dropdims(weights(A); dims = 3)) == @inferred(weights(B))
        @test @inferred(dropdims(nobs(A); dims = 3)) == @inferred(nobs(B))
        @test @inferred(dropdims(mean(A); dims = 3)) == mean(B)

        A = IndependentStatistic(Float64, 1, Float64, size(A))
        @inferred fit!(A, Float32.(x), w)
        @test nobs(A) ≈ sum(w, dims = dims)

        A = IndependentStatistic(Float64, 1, Float64, size(A))
        @inferred fit!(A, x, 1)

        B = IndependentStatistic(Float64, 1, size(A))
        @inferred fit!(B, x, 1)

        C = IndependentStatistic(Float64, 1, Int, size(A))
        @inferred fit!(C, x, trues(size(x)...))

        @test @inferred(weights(A)) == @inferred(weights(B))
        @test @inferred(nobs(A)) == @inferred(nobs(B))
        @test @inferred(mean(A)) == mean(B)
        @test @inferred(weights(A)) == @inferred(weights(C))
        @test @inferred(nobs(A)) == @inferred(nobs(C))
        @test @inferred(mean(A)) == mean(C)
    end

    @testset "fit! from IndependentStatistic" begin
        x = randn(2, 3, 20)
        x1 = x[:, :, 1:10]
        x2 = x[:, :, 11:20]

        A = IndependentStatistic(4, x1; dims = 3)
        B = IndependentStatistic(4, x2; dims = 3)
        C = IndependentStatistic(4, x; dims = 3)
        fit!(A, B)

        @test isapprox(mean(A), mean(C); rtol = 1.0e-10)
        @test isapprox(var(A; corrected = false), var(C; corrected = false); rtol = 1.0e-10)

        w = rand(2, 3, 20)
        w1 = w[:, :, 1:10]
        w2 = w[:, :, 11:20]

        Aw = IndependentStatistic(2, x1, w1; dims = 3)
        Bw = IndependentStatistic(2, x2, w2; dims = 3)
        Cw = IndependentStatistic(2, x, w; dims = 3)
        fit!(Aw, Bw)

        @test isapprox(mean(Aw), mean(Cw); rtol = 1.0e-10)
        @test isapprox(var(Aw; corrected = false), var(Cw; corrected = false); rtol = 1.0e-10)
    end

    @testset "merge and merge!" begin
        x = randn(2, 3, 20)
        A = IndependentStatistic(2, x[:, :, 1:10]; dims = 3)
        B = IndependentStatistic(2, x[:, :, 11:20]; dims = 3)
        C = IndependentStatistic(2, x; dims = 3)

        M = merge(A, B)
        @test isapprox(mean(M), mean(C); rtol = 1.0e-10)
        @test isapprox(var(M; corrected = false), var(C; corrected = false); rtol = 1.0e-10)

        A2 = IndependentStatistic(2, x[:, :, 1:10]; dims = 3)
        merge!(A2, B)
        @test isapprox(mean(A2), mean(C); rtol = 1.0e-10)
        @test isapprox(var(A2; corrected = false), var(C; corrected = false); rtol = 1.0e-10)

        w = rand(2, 3, 20)
        Aw = IndependentStatistic(2, x[:, :, 1:10], w[:, :, 1:10]; dims = 3)
        Bw = IndependentStatistic(2, x[:, :, 11:20], w[:, :, 11:20]; dims = 3)
        Cw = IndependentStatistic(2, x, w; dims = 3)

        Mw = merge(Aw, Bw)
        @test isapprox(mean(Mw), mean(Cw); rtol = 1.0e-10)
        @test isapprox(var(Mw; corrected = false), var(Cw; corrected = false); rtol = 1.0e-10)

        Aw2 = IndependentStatistic(2, x[:, :, 1:10], w[:, :, 1:10]; dims = 3)
        merge!(Aw2, Bw)
        @test isapprox(mean(Aw2), mean(Cw); rtol = 1.0e-10)
        @test isapprox(var(Aw2; corrected = false), var(Cw; corrected = false); rtol = 1.0e-10)
    end
end
