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
    @test all(hdu -> hdu.data_size == dims, moments_hdus)
    @test weights_hdu.data_size == dims
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "Base.write(moments_hdus, weights_hdu, stat)" begin
    data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    K = 5
    dims = (3, 3)
    stat = IndependentStatistic(K, dims)
    fit!(stat, data)
    f = openfits(joinpath(dir, "test.fits"), "w!")
    (moments_hdus, weights_hdu) = OnlineSampleStatistics.declare_HDUs(f, stat)
    group_id = 1
    @test_nowarn Base.write(moments_hdus, weights_hdu, stat; group_id)
    @test OnlineSampleStatistics.isa_stat_hdu(weights_hdu)
    @test all(hdu -> OnlineSampleStatistics.isa_stat_hdu(hdu), moments_hdus)
    @test weights_hdu[ext.STAT_WEIGHTS_KWD].logical
    @test weights_hdu[ext.STAT_GROUP_ID_KWD].integer == group_id
    @test all(hdu -> hdu[ext.STAT_GROUP_ID_KWD].integer == group_id, moments_hdus)
    @test read(weights_hdu) == fill(1, dims)
    @test read(moments_hdus[1]) == data
    @test read(moments_hdus[2]) == fill(0, dims)
    @test read(moments_hdus[K]) == fill(0, dims)
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "Base.write(fitsfile, header, data)" begin
    data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    K = 5
    dims = (3, 3)
    stat = IndependentStatistic(K, dims)
    fit!(stat, data)
    f = openfits(joinpath(dir, "test.fits"), "w!")
    hdr = FitsHeader("TESTKWD" => true)
    otherhdr = FitsHeader("OTHER" => 42)
    otherdata = fill(1, 10, 10)
    @test_nowarn Base.write(f, hdr, stat, otherhdr, otherdata)
    @test length(f) == K+2
    @test all(i -> f[i]["TESTKWD"].logical, 1:(K+1))
    @test f[K+2]["OTHER"].integer == 42
    @test read(f[K+2]) == otherdata
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "isa_stat_hdu" begin
    data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    K = 5
    dims = (3, 3)
    stat = IndependentStatistic(K, dims)
    fit!(stat, data)
    f = openfits(joinpath(dir, "test.fits"), "w!")
    hdr = FitsHeader("TESTKWD" => true)
    Base.write(f, hdr, stat)
    @test all(i -> OnlineSampleStatistics.isa_stat_hdu(f[i]), 1:K+1)
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "Base.read(IndependentStatistic, fitsfile, hdu_ext)" begin
    data = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    K = 5
    dims = (3, 3)
    stat = IndependentStatistic(K, dims)
    fit!(stat, data)
    f = openfits(joinpath(dir, "test.fits"), "w!")
    hdr = FitsHeader("TESTKWD" => true)
    Base.write(f, hdr, stat)
    close(f)
    f = openfits(joinpath(dir, "test.fits"), "r")
    local stat2
    @test_nowarn stat2 = Base.read(IndependentStatistic, f, 1)
    try close(f) catch err; end # finally close file ignoring any error
end

end
end

