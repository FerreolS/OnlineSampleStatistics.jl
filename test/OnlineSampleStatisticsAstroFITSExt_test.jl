using Test, OnlineSampleStatistics, AstroFITS

Ext = Base.get_extension(OnlineSampleStatistics, :OnlineSampleStatisticsAstroFITSExt)

@testset "OnlineSampleStatisticsAstroFITSExt" begin
mktempdir() do dir

fitsfile = openfits(joinpath(dir, "test.fits"), "w!")

T1 = Float32
T2 = Float64

N = 2
K = 5
dims = (3, 3)

stat1 = IndependentStatistic(T1, K, dims)
stat2 = IndependentStatistic(T2, K, dims)

fit!(stat1, rand(T1, dims...))
fit!(stat1, rand(T1, dims...))

fit!(stat2, rand(T2, dims...))
fit!(stat2, rand(T2, dims...))
fit!(stat2, rand(T2, dims...))

stat_group_id1 = "ABC"
stat_group_id2 = "DEF"

hdr = FitsHeader("COUNT" => 42)

@testset "write" begin
    @test_nowarn write(fitsfile, hdr, stat1, stat_group_id1)
    @test length(fitsfile) == K+1
    @test_nowarn write(fitsfile, hdr, stat2, stat_group_id2)
    @test length(fitsfile) == 2(K+1)
    @test fitsfile[1]["COUNT"].integer == 42
    @test !haskey(fitsfile[2], "COUNT")
    @test fitsfile[K+2]["COUNT"].integer == 42
end

@testset "find_stat_groupd_ids" begin
    @test_nowarn Ext.find_stat_groupd_ids(fitsfile)
    @test Ext.find_stat_groupd_ids(fitsfile) == [stat_group_id1, stat_group_id2]
end

@testset "isa_stat_hdu" begin
    @test_nowarn Ext.isa_stat_hdu(fitsfile[1])
    @test Ext.isa_stat_hdu(fitsfile[1])
end

@testset "find_stat_hdus" begin
    @test_nowarn Ext.find_stat_hdus(fitsfile, stat_group_id1)
    @test length(Ext.find_stat_hdus(fitsfile, stat_group_id1)[1]) == K
    @test all(Ext.isa_stat_hdu, Ext.find_stat_hdus(fitsfile, stat_group_id1)[1])
    @test Ext.isa_stat_hdu(Ext.find_stat_hdus(fitsfile, stat_group_id1)[2])
    @test Ext.find_stat_hdus(fitsfile, stat_group_id1)[3] == T1
    @test Ext.find_stat_hdus(fitsfile, stat_group_id1)[4] == N
    @test Ext.find_stat_hdus(fitsfile, stat_group_id1)[5] == K
    @test Ext.find_stat_hdus(fitsfile, stat_group_id1)[6] <: Integer
    @test_nowarn Ext.find_stat_hdus(fitsfile, stat_group_id2)
end

@testset "read" begin
    @test_nowarn read(IndependentStatistic, fitsfile, stat_group_id1)
    @test nobs(read(IndependentStatistic, fitsfile, stat_group_id1)) == nobs(stat1)
    @test all(1:K) do i
              get_moments(read(IndependentStatistic, fitsfile, stat_group_id1), i) == 
              get_moments(stat1, i)
          end
    @test_nowarn read(IndependentStatistic, fitsfile, stat_group_id2)
    @test nobs(read(IndependentStatistic, fitsfile, stat_group_id2)) == nobs(stat2)
    @test all(1:K) do i
              get_moments(read(IndependentStatistic, fitsfile, stat_group_id2), i) == 
              get_moments(stat2, i)
          end
end

try close(fitsfile) catch err; end # finally close file ignoring any error

end
end

