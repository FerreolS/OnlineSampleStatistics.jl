using ZippedArrays, StructuredArrays

@testset "IndependentStatistic Tests" begin
    # Test constructor with default type
    @testset "Constructor Tests" begin
        A = IndependentStatistic((3, 3), 5)
        @test size(A) == (3, 3)
        @test typeof(A) <: ZippedArray
    end

    # Test constructor with custom type
    @testset "other Type Constructor" begin
        A = IndependentStatistic(Float32, (4, 4), 2)
        @test size(A) == (4, 4)
        @test typeof(A) <: ZippedArray
    end

    # Test mean calculation
    @testset "Mean & variance Calculation" begin
        using OnlineSampleStatistics: get_rawmoments, weights, get_moments
        data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        A = IndependentStatistic((3, 3), 5)
        push!(A, data)
        @test A == IndependentStatistic(data, 5)
        @test @inferred(order(A)) == 5
        B = [UnivariateStatistic(d, 5) for d ∈ data]
        @test @inferred(order(B)) == 5
        @test A == B
        @test @inferred(get_rawmoments(A, 1)) == @inferred(get_rawmoments(B, 1))
        @test @inferred(get_rawmoments(A)) == @inferred(get_rawmoments(B))
        @test @inferred(get_moments(A, 1)) == @inferred(get_moments(B, 1))
        @test @inferred(weights(A)) == @inferred(weights(B))
        @test @inferred(nobs(A)) == @inferred(nobs(B))
        @test @inferred(mean(A)) == data

        @inferred(push!(A, -1 .* data))
        @test @inferred(var(A; corrected=false)) == data .^ 2
        @test @inferred(var(A)) == 2 .* data .^ 2
        @test @inferred(skewness(A)) == zeros(Float64, 3, 3)
        @test @inferred(kurtosis(A)) == -2.0 * ones(3, 3)
        @test @inferred(var(A)) ≈ var.(A)
        @test @inferred(kurtosis(A)) ≈ kurtosis.(A)
        @test @inferred(skewness(A)) ≈ skewness.(A)
        @test @inferred(mean(A)) ≈ mean.(A)
        @test @inferred(var(A; corrected=false)) ≈ var.(A; corrected=false)
    end


    # Test error handling for invalid dimensions
    @testset "Error Handling" begin
        A = IndependentStatistic((3, 3), 2)
        @test_throws DimensionMismatch push!(A, [1.0 2.0; 3.0 4.0])
    end

    @testset "Weighted Data" begin
        x = randn(2, 3, 10)
        A = IndependentStatistic(x, trues(size(x)...), 2)
        B = IndependentStatistic(x, 2)
        @test @inferred(weights(A)) == @inferred(weights(B))
        @test @inferred(nobs(A)) == @inferred(nobs(B))
        @test @inferred(mean(A)) == mean(B)


        dims = 3
        A = IndependentStatistic(x, trues(size(x)...), 2; dims=dims)
        B = IndependentStatistic(x, 2; dims=dims)
        @test @inferred(weights(A)) == @inferred(weights(B))
        @test @inferred(nobs(A)) == @inferred(nobs(B))
        @test @inferred(mean(A)) == mean(B)
        @test @inferred(var(A, corrected=false)) == var(B, corrected=false)


        w = (rand(size(x)...) .> 0.1)
        A = IndependentStatistic(x, w, 2; dims=dims)

        @test mean(A) ≈ sum(w .* x; dims=dims) ./ sum(w; dims=dims)
        @test var(A; corrected=false) ≈ sum(w .* (x .- mean(A)) .^ 2; dims=dims) ./ sum(w; dims=dims)

        A = IndependentStatistic(Float64, size(A), Float64, 1)
        @inferred push!(A, x, w)
        @test nobs(A) ≈ sum(w, dims=dims)
    end
end