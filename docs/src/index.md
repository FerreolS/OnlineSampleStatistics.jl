```@meta
CurrentModule = OnlineSampleStatistics
```

# OnlineSampleStatistics

Documentation for [OnlineSampleStatistics](https://github.com/FerreolS/OnlineSampleStatistics.jl).

```@index
```

```@autodocs
Modules = [OnlineSampleStatistics, OnlineSampleStatisticsAstroFITSExt]
```

## I/O with AstroFITS

OnlineSampleStatistics provides FITS input/output for `IndependentStatistic`
through an extension loaded when `AstroFITS` is available.

```julia
using OnlineSampleStatistics
using AstroFITS
```

Write an `IndependentStatistic` to a FITS file:

```julia
stat = IndependentStatistic(2, (3, 3))
fit!(stat, rand(3, 3))

fitsfile = openfits("myfile.fits", "w!")
write(fitsfile, FitsHeader(), stat)
close(fitsfile)
```

Read it back:

```julia
fitsfile = openfits("myfile.fits")
stat = read(IndependentStatistic, fitsfile)
close(fitsfile)
```

Each raw moment is stored in its own HDU, and weights (`nobs`) are stored in a
dedicated HDU.

If you store multiple statistics in one file, use group IDs:

```julia
fitsfile = openfits("myfile.fits", "w!")

stat1 = IndependentStatistic(2, (3, 3))
fit!(stat1, rand(3, 3))
write(fitsfile, FitsHeader(), stat1, "ABC")

stat2 = IndependentStatistic(2, (3, 3))
fit!(stat2, rand(3, 3))
write(fitsfile, FitsHeader(), stat2, "DEF")

close(fitsfile)
```

```julia
fitsfile = openfits("myfile.fits")
stat1 = read(IndependentStatistic, fitsfile, "ABC")
stat2 = read(IndependentStatistic, fitsfile, "DEF")
close(fitsfile)
```

If no group ID is provided, one is generated automatically and checked for
uniqueness in the file. If you provide a group ID that already exists, writing
throws an `ArgumentError`.

To inspect available group IDs:

```julia
fitsfile = openfits("myfile.fits")
ids = OnlineSampleStatistics.find_stat_group_ids(fitsfile)
close(fitsfile)
```
