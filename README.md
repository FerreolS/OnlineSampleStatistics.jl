# OnlineSampleStatistics

[![Documentation][doc-img]][doc-url][![License][license-img]][license-url] [![Build Status][github-ci-img]][github-ci-url] [![Coverage][codecov-img]][codecov-url] [![Aqua QA][aqua-img]][aqua-url]

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[github-ci-img]: https://github.com/FerreolS/OnlineSampleStatistics.jl/actions/workflows/CI.yml/badge.svg?branch=master
[github-ci-url]: https://github.com/FerreolS/OnlineSampleStatistics.jl/actions/workflows/CI.yml?query=branch%3Amaster
[codecov-img]: http://codecov.io/github/FerreolS/OnlineSampleStatistics.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/FerreolS/OnlineSampleStatistics.jl?branch=master
[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl
[doc-img]: https://img.shields.io/badge/docs-latest-blue.svg
[doc-url]: https://ferreols.github.io/OnlineSampleStatistics.jl/dev/

`OnlineSampleStatistics.jl` is a Julia package for online, single-pass estimation of statistical moments. It implements formulae from [Pébay et al 2016](https://doi.org/10.1007/s00180-015-0637-z)

## Features

- compute any moment (mean, variance, skewness, kurtosis, etc.) in a single pass
- numerically stable (avoids catastrophic cancellation even for non-centered data)
- memory efficient (0 allocations) and fast ($\approx 10\\,\text{ns}$  per sample for the first two moments)
- handle weighted data
- cope with univariate and multivariate (array) data
  
Designed for scenarios where data arrives in a streaming fashion or when memory efficiency (0 allocations) is critical.

## Usage

It mainly provides two types: `UnivariateStatistic` and `IndependentStatistic`. For both types, samples are added using `fit!` methods, including weighted samples. `Statistics` and `StatsBase` methods are used to query `mean`, `var`, `std`, `skewness`, `kurtosis`, `wsum`, `weights`, and `nobs`.

`UnivariateStatistic` and `IndependentStatistic` objects support pretty-printing via `show` methods.


### Univariate Statistics

`UnivariateStatistic{T,K}` tracks `K` statistical moments (of type `T`) of a data stream. It is defined as a subtype of `OnlineStat{T}` from [OnlineStats.jl](https://github.com/joshday/OnlineStats.jl) to leverage its functionality, including [Transducers.jl](https://github.com/JuliaFolds/Transducers.jl) methods for parallel processing.

Constructor summary (preferred API):
- `UnivariateStatistic(K)`
- `UnivariateStatistic(K, x)`
- `UnivariateStatistic(K, x, w)`
- `UnivariateStatistic(K, x_array, w_array)`
- `UnivariateStatistic(T, K)` / `UnivariateStatistic(T, K, I)`
- `UnivariateStatistic(T, K, I, x, w)` for fully explicit typing

where `K` is the number of tracked moments, `I` is the weight type, `x` is a sample,
and `w` a weight. In `K`-first constructors, value type is inferred from `x`
(integer inputs are promoted to `Float64`).

```julia-repl
julia> using OnlineSampleStatistics

julia> # Create a univariate statistic to track up to 4 moments
julia> stat = UnivariateStatistic(4)
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 0


julia> # Add data incrementally
julia> fit!(stat, 1.0)                                             # single sample
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 1
  μ: 1.0
  σ²: 0.0  σ: 0.0
  skewness: NaN
  kurtosis: NaN

julia> fit!(stat, [2.0, 3.0, 4.0])                                # array of samples
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 3
  μ: 3.0
  σ²: 0.67  σ: 0.82
  skewness: 0.0
  kurtosis: -1.5

julia> fit!(stat, skipmissing( [missing, 5.0, 6.0, 7.0, missing])) # iterator of samples
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 6
  μ: 4.5
  σ²: 2.9  σ: 1.7
  skewness: -1.2e-16
  kurtosis: -1.3

julia> # Compute statistics
julia> println("Number of samples: ", nobs(stat))
Number of samples: 7

julia> println("Mean: ", mean(stat))
Mean: 4.0

julia> println("Variance: ", var(stat))
Variance: 4.666666666666667

julia> println("Skewness: ", skewness(stat))
Skewness: -1.2688263138573217e-16

julia> println("Kurtosis: ", kurtosis(stat))
Kurtosis: -1.2500000000000002

julia> # Compute weighted statistics
julia> fit!(stat, 1.0,2.0)                           # single sample with weight
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 9
  μ: 3.3
  σ²: 4.7  σ: 2.2
  skewness: 0.36
  kurtosis: -1.3

julia> fit!(stat, [2.0, 3.0, 4.0], [1.0, 2.0, 3.0])  # array of samples with weights
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 15
  μ: 3.3
  σ²: 3.0  σ: 1.7
  skewness: 0.39
  kurtosis: -0.55
```

### Independent Statistics

`IndependentStatistic{T,N,K}` tracks `K` statistical moments of an independent multivariate data stream. It is a subtype of `AbstractArray{UnivariateStatistic{T,K},N}` that uses a `ZippedArrays` to ensure the efficiency of the updates.

```julia-repl
julia> using OnlineSampleStatistics

julia> x = randn(3, 4, 10); #  array of samples

julia> w = rand(3, 4, 10);  #  array of weights

julia> stat = IndependentStatistic(2, x, w; dims=3); # Tracking 2 statistical moments along dimension dims=3.

julia> summary(stat)
"3×4×1 IndependentStatistic{Float64, 3, 2, Float64} with 2 moments"

julia> wmean = mean(stat);           # weighted mean

julia> summary(wmean)
"3×4×1 Array{Float64, 3}"

```

## Comparison with [OnlineStats.jl](https://github.com/joshday/OnlineStats.jl)

### Precision

Precisions are similar for both packages on centered data. However, `OnlineStats.jl` is less  numerically stable and suffers from catastrophic precision issues on non-centered data (leading to negative variance).

#### Centered data

```julia-repl

julia> N = 100_000_000;
julia> x =  (randn(N)) .^ 2;

julia> A = UnivariateStatistic(4); # accumulating moments up to 4th
julia> fit!(A,x);

julia> m = Moments();
julia> fit!(m,x);

julia> mean(m) , mean(A), mean(x)
(1.0001326319698987, 1.0001326319698987, 1.0001326319702761)

julia> var(m) , var(A), var(x)
(1.9996036953234633, 1.999603695323393, 1.9996036953236203)

julia> skewness(m), skewness(A), skewness(x)
(2.825531520993099, 2.8255315209913205, 2.8255315209922784)

julia> kurtosis(m) , kurtosis(A), kurtosis(x)
(11.961958271406537, 11.961958260570487, 11.961958271313998)

```

#### Non centered data

```julia-repl
julia> N = 100_000_000;
julia> x = 1e9 .+ (randn(N)) .^ 2;

julia> A = UnivariateStatistic(4); # accumulating moments up to 4th
julia> fit!(A,x);

julia> m = Moments();
julia> fit!(m,x);

julia> mean(m) , mean(A), mean(x)
(1.0000000014725969e9, 1.0000000014725969e9, 1.0000000009999052e9)

julia> var(m) , var(A), var(x)
(-1.7655987376559874e8, 2.1162718296263807, 1.9991845743442014)

julia> skewness(m)
ERROR: DomainError with -1.76559872e8:
Exponentiation yielding a complex result requires a complex argument.
Replace x^y with (x+0im)^y, Complex(x)^y, or similar.
Stacktrace:
 [1] throw_exp_domainerror(x::Float64)
   @ Base.Math ./math.jl:41
 [2] ^
   @ ./math.jl:1157 [inlined]
 [3] skewness(o::Moments{EqualWeight})
   @ OnlineStatsBase ~/.julia/packages/OnlineStatsBase/pYmRb/src/stats.jl:522
 [4] top-level scope
   @ REPL[121]:1

julia> skewness(A), skewness(x)
(2.599982413924975, 2.827216059756395)

julia> kurtosis(m) , kurtosis(A), kurtosis(x)
(1.5363254476149368e9, 10.688829742853885, 11.981695139465826)

```

### Performance

Performance is similar for both packages for the first two moments. However the additional computation in `OnlineSampleStatistics.jl` to ensure numerical precision for higher moments scales exponentially with the number of stored moments.

#### Mean

```julia-repl
julia> N = 100_000_000;
julia> x =  (randn(N)) .^ 2;
julia> m = Mean();
julia> M = UnivariateStatistic(1); # accumulating mean only (first moment)

julia> @btime fit!(m,x)
  274.430 ms (0 allocations: 0 bytes)

julia> @btime fit!(M,x);
  280.195 ms (0 allocations: 0 bytes)
```

#### Variance

```julia-repl
julia> v = Variance();
julia> A = UnivariateStatistic(2); 

julia> @btime fit!(v,x);
  311.701 ms (0 allocations: 0 bytes)

julia> @btime fit!(A,x);
  288.229 ms (0 allocations: 0 bytes)
```

#### Higher moments

```julia-repl
julia> m = Moments();

julia> @btime fit!(m,x)
  504.927 ms (0 allocations: 0 bytes)

julia> A = UnivariateStatistic(3); # accumulating moments up to 3th
julia> @btime fit!(A,x)
  575.805 ms (0 allocations: 0 bytes)

julia> A = UnivariateStatistic(4); # accumulating moments up to 4th
julia> @btime fit!(A,x)
  2.135 s (0 allocations: 0 bytes)

julia> A = UnivariateStatistic(5); # accumulating moments up to 5th
julia> @btime fit!(A,x)
  3.898 s (0 allocations: 0 bytes)
```

### Input Output with FITS files

Need to load package [AstroFITS](https://github.com/emmt/AstroFITS.jl) so the extension is loaded too:
```
using AstroFITS
```

Write an `IndependentStatistic` to a FITS file:
```
stat = IndependentStatistic(2, (3, 3))
fit!(stat, rand(3, 3))

fitsfile = openfits("myfile.fits", "w!")
write(fitsfile, FitsHeader(), stat)
close(fitsfile)
```

Read an `IndependantStatistic` from a FITS file:
```
fitsfile = openfits("myfile.fits")
stat = read(IndependentStatistic, fitsfile)
close(fitsfile)
```

Each (raw) moment is written in its own HDU. the weights (the `nobs`) are written in their own HDU too.

If you store multiple statistics inside the same FITS file, "group IDs" are used to identify which HDU belongs to which stat:

```
group_id1 = "ABC"
stat1 = IndependentStatistic(2, (3, 3))
fit!(stat1, rand(3, 3))

group_id2 = "DEF"
stat2 = IndependentStatistic(2, (3, 3))
fit!(stat2, rand(3, 3))

fitsfile = openfits("myfile.fits", "w!")
write(fitsfile, FitsHeader(), stat1, group_id1)
write(fitsfile, FitsHeader(), stat2, group_id2)
close(fitsfile)
```

```
fitsfile = openfits("myfile.fits")
stat1 = read(IndependentStatistic, fitsfile, "ABC")
stat2 = read(IndependentStatistic, fitsfile, "DEF")
close(fitsfile)
```


