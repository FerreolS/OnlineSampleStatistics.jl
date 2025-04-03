using ZippedArrays, StructuredArrays

@testset "IndependentStatistic Tests" begin
    # Test constructor with default type
    @testset "Constructor Tests" begin
        A = IndependentStatistic(2, (3, 3))
        @test size(A) == (3, 3)
        @test typeof(A) <: ZippedArray
    end

    # Test constructor with custom type
    @testset "Custom Type Constructor" begin
        A = IndependentStatistic(Float32, 2, (4, 4))
        @test size(A) == (4, 4)
        @test typeof(A) <: ZippedArray
    end

    # Test mean calculation
    @testset "Mean & variance Calculation" begin
        A = IndependentStatistic(2, (3, 3))
        push!(A, [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])
        @test mean(A) == [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

        @inferred(push!(A, -1 .* [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]))
        @test @inferred(var(A; corrected=false)) == [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0] .^ 2
    end

    #=     # Test skewness calculation
        @testset "Skewness Calculation" begin
            A = IndependentStatistic(3, (3, 3))
            push!(A, [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])
            @test StatsBase.skewness(A) == [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        end

        # Test kurtosis calculation
        @testset "Kurtosis Calculation" begin
            A = IndependentStatistic(4, (3, 3))
            push!(A, [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])
            @test StatsBase.kurtosis(A) == [-3.0 -3.0 -3.0; -3.0 -3.0 -3.0; -3.0 -3.0 -3.0]
        end =#

    # Test error handling for invalid dimensions
    @testset "Error Handling" begin
        A = IndependentStatistic(2, (3, 3))
        @test_throws DimensionMismatch push!(A, [1.0 2.0; 3.0 4.0])
    end
end