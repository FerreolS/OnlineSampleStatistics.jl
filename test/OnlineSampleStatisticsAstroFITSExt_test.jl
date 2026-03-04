using Test, OnlineSampleStatistics, AstroFITS, FITSHeaders

ext = Base.get_extension(OnlineSampleStatistics, :OnlineSampleStatisticsAstroFITSExt)

@testset "OnlineSampleStatisticsAstroFITSExt" begin
mktempdir() do dir

@testset "declare_HDUs" begin
    data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    K = 5
    dims = (3, 3)
    stat = IndependentStatistic(K, dims)
    fit!(stat, data)
    f = openfits(joinpath(dir, "test.fits"), "w!")
    local moments_hdus, weights_hdu
    @test_nowarn (moments_hdus, weights_hdu) = OnlineSampleStatistics.declare_HDUs(f, stat)
    @test length(moments_hdus) == K
    @test all(hdu -> size(hdu) == dims, moments_hdus)
    @test size(weights_hdu) == dims
    try close(f) catch err; end # finally close file ignoring any error
end

end
end
