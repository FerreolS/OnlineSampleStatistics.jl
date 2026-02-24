using StructuredArrays
using Statistics

"""
    show(io::IO, A::UnivariateStatistic)

Pretty-print a `UnivariateStatistic` object in a readable format with statistics summary.

Displays:
- Number of observations (n)
- Order (number of moments tracked)
- Mean and variance (if available)
- Skewness and kurtosis (if K ≥ 3 or K ≥ 4)
"""
function Base.show(io::IO, ::MIME"text/plain", A::UnivariateStatistic{T, K, I}) where {T, K, I}
    n = nobs(A)
    println(io, "UnivariateStatistic{$T, $K, $I} with $K moments")
    if I <: Integer
        println(io, "  nobs: $n")
    else
        println(io, "  weight: $(round(n; sigdigits = 3))")
    end

    if n == 0
        return
    end

    # Always show mean
    μ = mean(A)
    print(io, "  μ: ", sprint(show, μ))

    # Show variance and std if K ≥ 2
    if K ≥ 2
        σ² = var(A; corrected = false)
        σ = std(A; corrected = false)
        println(io)
        print(io, "  σ²: ", sprint(show, σ²), "  σ: ", sprint(show, σ))
    end

    # Show skewness if K ≥ 3
    if K ≥ 3
        γ = skewness(A)
        println(io)
        print(io, "  skewness: ", sprint(show, γ))
    end

    # Show kurtosis if K ≥ 4
    if K ≥ 4
        κ = kurtosis(A)
        println(io)
        print(io, "  kurtosis: ", sprint(show, κ))
    end
    return
end

# Compact display for REPL and inline use
Base.show(io::IO, A::UnivariateStatistic) = Base.summary(io, A)

# One-line summary for compact listings (e.g. varinfo)
"""
    summary(A::UnivariateStatistic)

Return a compact one-line summary of the statistic.
"""
function Base.summary(io::IO, A::UnivariateStatistic{T, K, I}) where {T, K, I}
    return print(io, "UnivariateStatistic{$T,$K,$I}  with $K moments")
end

# Pretty-print for IndependentStatistic
function Base.show(io::IO, ::MIME"text/plain", A::IndependentStatistic{T, N, K, I}) where {T, N, K, I}
    wts = weights(A)
    is_uniform = wts isa MutableUniformArray

    if is_uniform
        println(io, "$(Base.dims2string(size(A))) IndependentStatistic{$T, $N, $K, $I} with $K moments")
        println(io, "  nobs: $(StructuredArrays.value(wts))")
    else
        println(io, "$(Base.dims2string(size(A)))  Weighted IndependentStatistic{$T, $N, $K, $I} with $K moments")

        # Compute weight statistics with reduced precision
        wts_flat = vec(wts)
        wts_min = round(minimum(wts_flat); sigdigits = 3)
        wts_max = round(maximum(wts_flat); sigdigits = 3)
        wts_median = round(Statistics.median(wts_flat); sigdigits = 3)
        println(io, "weights: min=$wts_min, median=$wts_median, max=$wts_max")
    end
    return
end

"""
    summary(A::IndependentStatistic)

Return a compact one-line summary of the statistic.
"""
function Base.summary(io::IO, A::IndependentStatistic{T, N, K, I}) where {T, N, K, I}
    return print(io, "$(Base.dims2string(size(A))) IndependentStatistic{$T, $N, $K, $I} with $K moments")
end

# Compact display for REPL and inline use
Base.show(io::IO, A::IndependentStatistic) = Base.summary(io, A)
