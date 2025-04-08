using ZippedArrays, StructuredArrays

@testset "IndependentStatistic Tests" begin
    # Test constructor with default type
    @testset "Constructor Tests" begin
        A = IndependentStatistic(2, (3, 3))
        @test size(A) == (3, 3)
        @test typeof(A) <: ZippedArray
    end

    # Test constructor with custom type
    @testset "other Type Constructor" begin
        A = IndependentStatistic(Float32, 2, (4, 4))
        @test size(A) == (4, 4)
        @test typeof(A) <: ZippedArray
    end

    # Test mean calculation
    @testset "Mean & variance Calculation" begin
        using OnlineSampleStatistics: get_rawmoments, weights, get_moments
        data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        A = IndependentStatistic(5, (3, 3))
        push!(A, data)
        @test A == IndependentStatistic(5, data)
        B = [UnivariateStatistic(d, 5) for d ∈ data]
        @test A == B
        @test get_rawmoments(A, 1) == get_rawmoments(B, 1)
        @test weights(A) == weights(B)
        @test nobs(A) == nobs(B)
        @test mean(A) == data

        @inferred(push!(A, -1 .* data))
        @test @inferred(var(A; corrected=false)) == data .^ 2
        @test @inferred(var(A)) == 2 .* data .^ 2
        @test @inferred(skewness(A)) == zeros(Float64, 3, 3)
        @test @inferred(kurtosis(A)) == -2.0 * ones(3, 3)
        @test var(A) ≈ var.(A)
        @test kurtosis(A) ≈ kurtosis.(A)
        @test skewness(A) ≈ skewness.(A)
        @test mean(A) ≈ mean.(A)
        @test var(A; corrected=false) ≈ var.(A; corrected=false)
    end


    # Test error handling for invalid dimensions
    @testset "Error Handling" begin
        A = IndependentStatistic(2, (3, 3))
        @test_throws DimensionMismatch push!(A, [1.0 2.0; 3.0 4.0])
    end

    @testset "Weighted Data" begin
        dims = 3
        x = randn(2, 3, 10)
        w = (rand(size(x)...) .> 0.1)
        A = IndependentStatistic(2, x, w, dims=dims)

        @test mean(A) ≈ sum(w .* x; dims=dims) ./ sum(w; dims=dims)
        @test var(A; corrected=false) ≈ sum(w .* (x .- mean(A)) .^ 2; dims=dims) ./ sum(w; dims=dims)

        A = IndependentStatistic(Float64, 1, size(A), Float64)
        @inferred push!(A, x, w)
        @test nobs(A) ≈ sum(w, dims=dims)
    end
end