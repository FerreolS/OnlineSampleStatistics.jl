# Core Types

## UnivariateStatistic

`UnivariateStatistic{T,K}` tracks online moments of a scalar stream of type `T` up to order `K`. 
Common constructors:
- `UnivariateStatistic(T, K, W)` empty statistic with value type `T`, moment order `K`, and weight type `W`;
- `UnivariateStatistic(T, K, W, x, w)` initialize from one weighted sample with types specified converting the sample if needed;
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

julia> cpx = 2. *cis.(π/3:π/3:2π)
6-element Vector{ComplexF64}:
  1.0000000000000002 + 1.7320508075688772im
 -0.9999999999999996 + 1.7320508075688774im
                -2.0 + 2.4492935982947064e-16im
 -1.0000000000000009 - 1.732050807568877im
  0.9999999999999987 - 1.732050807568878im
                 2.0 - 2.266215559059192e-15im

julia> cpxstat =UnivariateStatistic(2,cpx)
UnivariateStatistic{ComplexF64, 2, Int64} with 2 moments
  nobs: 6
  μ: -2.2e-16 - 4.2e-16im
  σ²: 4.0 + 0.0im  σ: 2.0 + 0.0im

julia> cpxstat =UnivariateStatistic(3,cpx)
ERROR: ArgumentError: UnivariateStatistic : 3 > 2 not implemented for complex numbers
```

### Weighted UnivariateStatistic

UnivariateStatistic supports weighted updates are supported with non-negative real weights.

Integer weights are usefull to discard some samples or to represents a number of observation.

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> weightedstat = UnivariateStatistic(2)
UnivariateStatistic{Float64, 2, Int64} with 2 moments
  nobs: 0

julia> fit!(weightedstat, [1.0, 3.0, -10.], [true, true, false])
UnivariateStatistic{Float64, 2, Int64} with 2 moments
  nobs: 2
  μ: 2.0
  σ²: 1.0  σ: 1.0

julia> fit!(weightedstat, [0, π, 2.18], [4, 1, 2])
UnivariateStatistic{Float64, 2, Int64} with 2 moments
  nobs: 9
  μ: 1.3
  σ²: 1.6  σ: 1.3

julia> mean(weightedstat)
1.2779547392877548

julia> var(weightedstat)
1.8344861950096323

julia> var(weightedstat; corrected=false)
1.6306543955641175
```

Real weights can represent sampling probability for each observation.

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> weightedstat = UnivariateStatistic(2,  [0, π, 2.18], [4.0, 1.0, 2.0])
UnivariateStatistic{Float64, 2, Float64} with 2 moments
  weight: 7.0
  μ: 1.1
  σ²: 1.6  σ: 1.3

julia> fit!(weightedstat, [0, 1, 2.1], [0.1, 0.1, 0.2])
UnivariateStatistic{Float64, 2, Float64} with 2 moments
  weight: 7.4
  μ: 1.1
  σ²: 1.6  σ: 1.3

julia> mean(weightedstat)
1.083999007241864

julia> var(weightedstat; corrected=false)
1.5758116119053236
  
```
!!! warning
  For consistency with other packages, `var` applies Bessel's correction when weights are integers (unbiased estimator) but not when weights are floats (frequency weights). Use the `corrected` keyword argument to precise this behavior.


## IndependentStatistic

`IndependentStatistic{T,N,K,W}` is an array-like container of independent univariate statistics.
It is useful for data where moments are accumulated independently per position (e.g. data from array of sensors and camera)

Typical workflow:
- build from an `N`-D sample array;
- optionally reduce over one or more sample dimensions via `dims`;
- query array-shaped `mean`, `var`, etc.

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> using OnlineSampleStatistics

julia> x = reshape(collect(1.0:120.0), 2, 3, 20);

julia> a = IndependentStatistic(2, x; dims=3);

julia> size(a), size(mean(a)), size(var(a))
((2, 3, 1), (2, 3, 1), (2, 3, 1))
```

Weighted construction is also supported:

```jldoctest; setup = :(using OnlineSampleStatistics)
julia> x = reshape(collect(1.0:120.0), 2, 3, 20);

julia> w = fill(2.0, 2, 3, 20);

julia> aw = IndependentStatistic(2, x, w; dims=3);

julia> size(mean(aw)), size(var(aw))
((2, 3, 1), (2, 3, 1))

julia> all(isfinite, var(aw))
true
```
