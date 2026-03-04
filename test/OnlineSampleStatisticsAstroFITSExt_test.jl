using Test, OnlineSampleStatistics, AstroFITS

astroext = Base.get_extension(OnlineSampleStatistics, :OnlineSampleStatisticsAstroFITSExt)

@testset "OnlineSampleStatisticsAstroFITSExt" begin
mktempdir() do dir

K = 5
dims = (3, 3)

stat1 = IndependentStatistic(K, dims)
stat2 = IndependentStatistic(K, dims)

fit!(stat1, rand(dims...))
fit!(stat1, rand(dims...))

fit!(stat2, rand(dims...))
fit!(stat2, rand(dims...))

group_id1 = "ABC"
group_id2 = "DEF"

@testset "declare_HDUs" begin
    f = openfits(joinpath(dir, "test.fits"), "w!")
    local moments_hdus, weights_hdu
    @test_nowarn (moments_hdus, weights_hdu) = astroext.declare_HDUs(f, stat1)
    @test length(moments_hdus) == K
    @test all(hdu -> hdu.data_size == dims, moments_hdus)
    @test weights_hdu.data_size == dims
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "Base.write(moments_hdus, weights_hdu, stat, group_id)" begin
    f = openfits(joinpath(dir, "test.fits"), "w!")
    (moments_hdus, weights_hdu) = astroext.declare_HDUs(f, stat1)
    @test_nowarn Base.write(moments_hdus, weights_hdu, stat1, group_id1)
    @test astroext.isa_stat_hdu(weights_hdu)
    @test all(hdu -> astroext.isa_stat_hdu(hdu), moments_hdus)
    @test weights_hdu[astroext.STAT_WEIGHTS_KWD].logical
    @test weights_hdu[astroext.STAT_GROUP_ID_KWD].string == group_id1
    @test all(hdu -> hdu[astroext.STAT_GROUP_ID_KWD].string == group_id1, moments_hdus)
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "Base.write(fitsfile, header, data, ?group_id)" begin
    f = openfits(joinpath(dir, "test.fits"), "w!")
    hdr = FitsHeader("TESTKWD" => true)
    @test_nowarn Base.write(f, hdr, stat1, group_id1)
    @test length(f) == K+1
    @test f[1]["TESTKWD"].logical
    @test f[1][astroext.STAT_GROUP_ID_KWD].string == group_id1
    @test !haskey(f[2], "TESTKWD")
    @test_nowarn Base.write(f, hdr, stat1)
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "isa_stat_hdu" begin
    f = openfits(joinpath(dir, "test.fits"), "w!")
    hdr = FitsHeader("TESTKWD" => true)
    Base.write(f, hdr, stat1)
    @test all(i -> astroext.isa_stat_hdu(f[i]), 1:K+1)
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "find_stat_groupd_ids" begin
    f = openfits(joinpath(dir, "test.fits"), "w!")
    Base.write(f, (), stat1, group_id1)
    Base.write(f, (), stat2, group_id2)
    @test OnlineSampleStatistics.find_stat_groupd_ids(f) == Set((group_id1, group_id2))
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "Base.read(IndependentStatistic, fitsfile, group_id)" begin
    f = openfits(joinpath(dir, "test.fits"), "w!")
    hdr = FitsHeader("TESTKWD" => true)
    Base.write(f, hdr, stat1)
    group_id = f[1][astroext.STAT_GROUP_ID_KWD].string
    close(f)
    f = openfits(joinpath(dir, "test.fits"), "r")
    local stat2
    @test_nowarn stat2 = Base.read(IndependentStatistic, f, group_id)
    try close(f) catch err; end # finally close file ignoring any error
end

@testset "Base.read(IndependentStatistic, fitsfile; ext)" begin
    f = openfits(joinpath(dir, "test.fits"), "w!")
    hdr = FitsHeader("TESTKWD" => true)
    Base.write(f, hdr, stat1)
    close(f)
    f = openfits(joinpath(dir, "test.fits"), "r")
    local stat1read
    @test_nowarn statread = Base.read(IndependentStatistic, f)
    @test all(i -> get_moments(stat1, i) == get_moments(statread, i), 1:K)
    @test nobs(stat1) == nobs(statread)
    @test_nowarn Base.read(IndependentStatistic, f; ext=1)
    @test_nowarn Base.read(IndependentStatistic, f; ext=2)
    try close(f) catch err; end # finally close file ignoring any error
end

end
end

