# Core Types

## UnivariateStatistic

`UnivariateStatistic{T,K}` tracks online moments of a scalar stream of type `T` up to order `K`.

Common constructors:
- `UnivariateStatistic(K)` empty statistic with `Float64` values and `Int` weights;
- `UnivariateStatistic(K, x)` initialize from one sample (type inferred from `x`);
- `UnivariateStatistic(K, x, w)` initialize from one weighted sample;
- `UnivariateStatistic(K, xs, ws)` initialize from weighted arrays.

Typical workflow:
- create with desired moment order `K`;
- update incrementally with `fit!`;
- query `nobs`, `mean`, `var`, `std`, `skewness`, `kurtosis` (when available).

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> using OnlineSampleStatistics

julia> s = UnivariateStatistic(4)
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 0

julia> fit!(s, [1.0, 2.0, 3.0, 4.0])
UnivariateStatistic{Float64, 4, Int64} with 4 moments
  nobs: 4
  μ: 2.5
  σ²: 1.2  σ: 1.1
  skewness: 0.0
  kurtosis: -1.4

julia> nobs(s)
4

julia> mean(s), var(s), skewness(s), kurtosis(s)
(2.5, 1.6666666666666667, 0.0, -1.36)

julia> s2 = UnivariateStatistic(Float64, 2, Float64, 1.5, 1.0)
UnivariateStatistic{Float64, 2, Float64} with 2 moments
  weight: 1.0
  μ: 1.5
  σ²: 0.0  σ: 0.0

julia> s3 = UnivariateStatistic(2, 1.5, 3.0)
UnivariateStatistic{Float64, 2, Float64} with 2 moments
  weight: 3.0
  μ: 1.5
  σ²: 0.0  σ: 0.0
```
### Weighted UnivariateStatistic

Weighted updates are supported. As an example these weights are usefull to represent the number of observation or sampling probability for each observation. In that case the moment weights are taken into an account in:

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> sw = UnivariateStatistic(2);

julia> fit!(sw, [1.0, 3.0, -10.], [1, 3, 0])
UnivariateStatistic{Float64, 2, Int64} with 2 moments
  nobs: 4
  μ: 2.5
  σ²: 0.75  σ: 0.87

julia> mean(sw)
2.5

julia> var(sw)
1.0

julia> var(sw; corrected=false)
0.75

```

## IndependentStatistic

`IndependentStatistic{T,N,K,W}` is an array-like container of independent univariate statistics.
It is useful for tensor data where moments are accumulated independently per position.

Typical workflow:
- build from an `N`-D sample array;
- optionally reduce over one or more sample dimensions via `dims`;
- query array-shaped `mean`, `var`, etc.

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> using OnlineSampleStatistics

julia> x = reshape(collect(1.0:120.0), 2, 3, 20);

julia> a = IndependentStatistic(x, 2; dims=3);

julia> size(a), size(mean(a)), size(var(a))
((2, 3, 1), (2, 3, 1), (2, 3, 1))
```

Weighted construction is also supported:

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> x = reshape(collect(1.0:120.0), 2, 3, 20);

julia> w = fill(2.0, 2, 3, 20);

julia> aw = IndependentStatistic(x, w, 2; dims=3);

julia> size(mean(aw)), size(var(aw))
((2, 3, 1), (2, 3, 1))

julia> all(isfinite, var(aw))
true
```
