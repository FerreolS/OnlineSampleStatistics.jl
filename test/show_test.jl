@testset "Show methods" begin
    
    @testset "UnivariateStatistic plain text show" begin
        A = UnivariateStatistic(2)
        fit!(A, [1.0, 2.0, 3.0])
        io = IOBuffer()
        show(io, MIME"text/plain"(), A)
        output = String(take!(io))

        @test contains(output, "UnivariateStatistic")
        @test contains(output, "nobs: 3")
        @test contains(output, "μ:")
        @test contains(output, "σ²:")
        @test contains(output, "σ:")
        @test contains(output, "with 2 moments")
    end

    @testset "UnivariateStatistic show with higher moments" begin
        A = UnivariateStatistic(4)
        fit!(A, randn(50))
        io = IOBuffer()
        show(io, MIME"text/plain"(), A)
        output = String(take!(io))

        @test contains(output, "skewness:")
        @test contains(output, "kurtosis:")
    end

    @testset "UnivariateStatistic empty show" begin
        A = UnivariateStatistic(3)
        io = IOBuffer()
        show(io, MIME"text/plain"(), A)
        output = String(take!(io))

        @test contains(output, "nobs: 0")
        @test !contains(output, "μ:")
    end

    @testset "UnivariateStatistic summary is one line" begin
        A = UnivariateStatistic(4)
        fit!(A, randn(10))

        io = IOBuffer()
        summary(io, A)
        output = String(take!(io))

        @test startswith(output, "UnivariateStatistic{")
        @test contains(output, "with 4 moments")
        @test !contains(output, "\n")
    end

    @testset "IndependentStatistic compact show" begin
        B = IndependentStatistic(randn(3, 3, 50), 4; dims = 3)
        io = IOBuffer()
        show(io, B)
        output = String(take!(io))

        @test contains(output, "IndependentStatistic")
        @test contains(output, "3×3×1")
    end

    @testset "IndependentStatistic summary is one line" begin
        B = IndependentStatistic(randn(3, 3, 20), 3; dims = 3)
        io = IOBuffer()
        summary(io, B)
        output = String(take!(io))

        @test contains(output, "3×3×1")
        @test contains(output, "with 3 moments")
        @test !contains(output, "\n")

        B_weighted = IndependentStatistic(randn(2, 4, 12), rand(2, 4, 12), 2; dims = 3)
        io = IOBuffer()
        summary(io, B_weighted)
        output = String(take!(io))

        @test contains(output, "2×4×1")
        @test contains(output, "with 2 moments")
        @test !contains(output, "\n")
    end

    @testset "IndependentStatistic plain text show" begin
        B = IndependentStatistic(randn(2, 3, 20), 2; dims = 3)
        io = IOBuffer()
        show(io, MIME"text/plain"(), B)
        output = String(take!(io))

        @test contains(output, "IndependentStatistic")
        @test contains(output, "2×3×1")
        @test contains(output, "with 2 moments")
    end

    @testset "Weighted IndependentStatistic show" begin
        B_weighted = IndependentStatistic(randn(2, 3, 20), rand(2, 3, 20), 2; dims = 3)
        io = IOBuffer()
        show(io, MIME"text/plain"(), B_weighted)
        output = String(take!(io))

        @test contains(output, "IndependentStatistic")
        @test contains(output, "2×3×1")
        @test contains(output, "with 2 moments")
    end
end
