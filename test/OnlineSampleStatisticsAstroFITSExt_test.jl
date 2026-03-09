using Test, OnlineSampleStatistics, AstroFITS, Random
using OnlineSampleStatistics: isa_stat_hdu, find_stat_group_ids

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
            @test length(fitsfile) == K + 1
            @test_nowarn write(fitsfile, hdr, stat2, stat_group_id2)
            @test length(fitsfile) == 2(K + 1)
            @test fitsfile[1]["COUNT"].integer == 42
            @test !haskey(fitsfile[2], "COUNT")
            @test fitsfile[K + 2]["COUNT"].integer == 42

            err = ArgumentError("stat group ID \"ABC\" already exists in FITS file")
            @test_throws err write(fitsfile, hdr, stat1, stat_group_id1)
        end

        @testset "write auto-generated group ID collision" begin
            collision_file = openfits(joinpath(dir, "collision.fits"), "w!")
            try
                # Build a deterministic candidate ID, write it once explicitly, then
                # reset RNG so the first auto-generated ID collides and must be retried.
                Random.seed!(1234)
                collision_id = string(rand('A':'Z', 16)...)

                @test_nowarn write(collision_file, FitsHeader(), stat1, collision_id)

                Random.seed!(1234)
                @test_nowarn write(collision_file, FitsHeader(), stat2)

                ids = find_stat_group_ids(collision_file)
                @test length(ids) == 2
                @test count(==(collision_id), ids) == 1
            finally
                close(collision_file)
            end
        end

        @testset "find_stat_group_ids" begin
            @test_nowarn find_stat_group_ids(fitsfile)
            @test find_stat_group_ids(fitsfile) == [stat_group_id1, stat_group_id2]
        end

        @testset "isa_stat_hdu" begin
            @test_nowarn isa_stat_hdu(fitsfile[1])
            @test isa_stat_hdu(fitsfile[1])
        end

        @testset "find_stat_hdus" begin
            @test_nowarn Ext.find_stat_hdus(fitsfile, stat_group_id1)
            @test length(Ext.find_stat_hdus(fitsfile, stat_group_id1)[1]) == K
            @test all(isa_stat_hdu, Ext.find_stat_hdus(fitsfile, stat_group_id1)[1])
            @test isa_stat_hdu(Ext.find_stat_hdus(fitsfile, stat_group_id1)[2])
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

            @test read(IndependentStatistic, fitsfile) == stat1
            @test read(IndependentStatistic, fitsfile; ext = 2) == stat1
            @test read(IndependentStatistic, fitsfile; ext = "WEIGHTS") == stat1
        end

        @testset "read errors" begin
            fitsfile2_path = joinpath(dir, "test2.fits")

            writefits!(fitsfile2_path, FitsHeader(), [0;;])
            err = ArgumentError("HDU \"1\" is not a stat HDU")
            @test_throws err FitsFile(f -> read(IndependentStatistic, f), fitsfile2_path)

            writefits!(
                fitsfile2_path,
                FitsHeader(),
                [0;;],
                filter(!is_structural, FitsHeader(fitsfile["MOMENT-1"])),
                read(fitsfile["MOMENT-1"])
            )
            err = ArgumentError("HDU \"1\" is not a stat HDU")
            @test_throws err FitsFile(f -> read(IndependentStatistic, f), fitsfile2_path)

            writefits!(
                fitsfile2_path,
                filter(!is_structural, FitsHeader(fitsfile["MOMENT-1"])),
                read(fitsfile["MOMENT-1"])
            )
            err = ArgumentError("could not find weights HDU")
            @test_throws err FitsFile(f -> read(IndependentStatistic, f), fitsfile2_path)

            writefits!(
                fitsfile2_path,
                filter(!is_structural, FitsHeader(fitsfile["WEIGHTS"])),
                read(fitsfile["WEIGHTS"])
            )
            err = ArgumentError("could not find any moment HDU")
            @test_throws err FitsFile(f -> read(IndependentStatistic, f), fitsfile2_path)

            writefits!(
                fitsfile2_path,
                filter(!is_structural, FitsHeader(fitsfile["WEIGHTS"])),
                read(fitsfile["WEIGHTS"]),
                filter(!is_structural, FitsHeader(fitsfile["MOMENT-1"])),
                read(fitsfile["MOMENT-1"]),
                filter(!is_structural, FitsHeader(fitsfile["MOMENT-1"])),
                read(fitsfile["MOMENT-1"])
            )
            err = ArgumentError("duplicate moment number 1 HDU for group \"ABC\"")
            @test_throws err FitsFile(f -> read(IndependentStatistic, f), fitsfile2_path)

            # malformed stat HDU with missing group id must fail with clear ArgumentError
            writefits!(
                fitsfile2_path,
                FitsHeader(
                    "STAT-HDU" => (true, "is a OnlineSampleStatistics.jl data"),
                    "STAT-NB-MOMENTS" => (K, "number of statistical moments"),
                    "STAT-MOMENT-INDEX" => (1, "th statistical moment")
                ),
                read(fitsfile["MOMENT-1"])
            )
            err = ArgumentError("HDU \"1\" is missing STAT-GROUP-ID keyword")
            @test_throws err FitsFile(f -> read(IndependentStatistic, f), fitsfile2_path)
        end

        try
            close(fitsfile)
        catch err
        end # finally close file ignoring any error

    end
end
