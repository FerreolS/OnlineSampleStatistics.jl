# OnlineSampleStatistics

[![License][license-img]][license-url] [![Build Status][github-ci-img]][github-ci-url] [![Coverage][codecov-img]][codecov-url] [![Aqua QA][aqua-img]][aqua-url]

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[github-ci-img]: https://github.com/FerreolS/OnlineSampleStatistics.jl/actions/workflows/CI.yml/badge.svg?branch=master
[github-ci-url]: https://github.com/FerreolS/OnlineSampleStatistics.jl/actions/workflows/CI.yml?query=branch%3Amaster
[codecov-img]: http://codecov.io/github/FerreolS/OnlineSampleStatistics.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/FerreolS/OnlineSampleStatistics.jl?branch=master
[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

## Comparison with [OnlineStats.jl](https://github.com/joshday/OnlineStats.jl)

### Precision

Precisions are similar for both packages on centered data. However, `OnlineStats.jl` is less  numerically stable and suffers from catastrophic precision issues on non-centered data (leading to negative variance).

#### Centered data

```julia-repl

julia> N = 100_000_000;
julia> x =  (randn(N)) .^ 2;

julia> A = UnivariateStatistic(4); # accumulating moments up to 4th
julia> push!(A,x);

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
julia> push!(A,x);

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
Performance is similar for both packages for the first two moments. However additional computations in `OninesSampleStatistics.jl` to ensure numerical precision for higher moments scale exponentially with the number of stored moment.


#### Mean

```julia-repl
julia> N = 100_000_000;
julia> x =  (randn(N)) .^ 2;
julia> m = Mean();
julia> M = UnivariateStatistic(1); # accumulating mean only (first moment)

julia> @btime fit!(m,x)
  274.430 ms (0 allocations: 0 bytes)

julia> @btime push!(A,x);
  280.195 ms (0 allocations: 0 bytes)
```

#### Variance

```julia-repl
julia> v = Variance();
julia> A = UnivariateStatistic(2); 

julia> @btime fit!(m,x)
  311.701 ms (0 allocations: 0 bytes)

julia> @btime push!(A,x);
  288.229 ms (0 allocations: 0 bytes)
```

#### Higher moments

```julia-repl
julia> m = Moments();
julia> A = UnivariateStatistic(3); # accumulating moments up to 3th
julia> @btime push!(A,x)
  575.805 ms (0 allocations: 0 bytes)

julia> @btime fit!(m,x)
  504.927 ms (0 allocations: 0 bytes)

julia> A = UnivariateStatistic(4); # accumulating moments up to 4th
julia> @btime push!(A,x)
  2.135 s (0 allocations: 0 bytes)

julia> A = UnivariateStatistic(5); # accumulating moments up to 5th
julia> @btime push!(A,x)
  3.898 s (0 allocations: 0 bytes)
```