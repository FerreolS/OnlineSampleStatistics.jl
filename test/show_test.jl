@testset "Show methods" begin
    @testset "UnivariateStatistic compact show" begin
        A = UnivariateStatistic(4)
        io = IOBuffer()
        show(io, A)
        output = String(take!(io))
        @test contains(output, "UnivariateStatistic")
        @test contains(output, "n=0")
        
        fit!(A, randn(100))
        io = IOBuffer()
        show(io, A)
        output = String(take!(io))
        @test contains(output, "UnivariateStatistic")
        @test contains(output, "100")
    end
    
    @testset "UnivariateStatistic plain text show" begin
        A = UnivariateStatistic(2)
        fit!(A, [1.0, 2.0, 3.0])
        io = IOBuffer()
        show(io, MIME"text/plain"(), A)
        output = String(take!(io))
        
        @test contains(output, "UnivariateStatistic")
        @test contains(output, "n: 3")
        @test contains(output, "μ:")
        @test contains(output, "σ²:")
        @test contains(output, "σ:")
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
        
        @test contains(output, "empty")
        @test !contains(output, "μ:")
    end
    
    @testset "IndependentStatistic compact show" begin
        B = IndependentStatistic(randn(3, 3, 50), 4; dims=3)
        io = IOBuffer()
        show(io, B)
        output = String(take!(io))
        
        @test contains(output, "IndependentStatistic")
        @test contains(output, "size=")
    end
    
    @testset "IndependentStatistic plain text show" begin
        B = IndependentStatistic(randn(2, 3, 20), 2; dims=3)
        io = IOBuffer()
        show(io, MIME"text/plain"(), B)
        output = String(take!(io))
        
        @test contains(output, "IndependentStatistic")
        @test contains(output, "size:")
        @test contains(output, "moments:")
        @test contains(output, "nobs:")
    end
    
    @testset "Weighted IndependentStatistic show" begin
        B_weighted = IndependentStatistic(randn(2, 3, 20), rand(2, 3, 20), 2; dims=3)
        io = IOBuffer()
        show(io, MIME"text/plain"(), B_weighted)
        output = String(take!(io))
        
        @test contains(output, "Weighted IndependentStatistic")
        @test contains(output, "size:")
        @test contains(output, "moments:")
        @test contains(output, "weights:")
        @test contains(output, "min=")
        @test contains(output, "median=")
        @test contains(output, "max=")
    end
end
