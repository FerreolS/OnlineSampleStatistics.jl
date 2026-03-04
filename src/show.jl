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
    sig_digits = get(io, :sig_digits, 2)

    n = nobs(A)
    println(io, "UnivariateStatistic{$T, $K, $I} with $K moments")
    if I <: Integer
        print(io, "  nobs: $n")
    else
        print(io, "  weight: $(round.(n; sigdigits = sig_digits))")
    end

    if n == 0
        return
    end

    # Always show mean
    μ = round(mean(A); sigdigits = sig_digits)
    println(io)
    print(io, "  μ: ", sprint(show, μ))

    # Show variance and std if K ≥ 2
    if K ≥ 2
        σ² = round(var(A; corrected = false); sigdigits = sig_digits)
        σ = round(std(A; corrected = false); sigdigits = sig_digits)
        println(io)
        print(io, "  σ²: ", sprint(show, σ²), "  σ: ", sprint(show, σ))
    end

    # Show skewness if K ≥ 3
    if K ≥ 3
        γ = round(skewness(A); sigdigits = sig_digits)
        println(io)
        print(io, "  skewness: ", sprint(show, γ))
    end

    # Show kurtosis if K ≥ 4
    if K ≥ 4
        κ = round(kurtosis(A); sigdigits = sig_digits)
        println(io)
        print(io, "  kurtosis: ", sprint(show, κ))
    end
    return
end

# Compact display for REPL and inline use
Base.show(io::IO, A::UnivariateStatistic{T, 1}) where {T} = print(io, mean(A))
function Base.show(io::IO, A::UnivariateStatistic{T, K}) where {T, K}
    error_digits = get(io, :error_digits, 2)
    return print(io, mean(A), " ± ", round(std(A; corrected = false), sigdigits = error_digits))
end


# One-line summary for compact listings (e.g. varinfo)
"""
	    summary(A::UnivariateStatistic)

Return a compact one-line summary of the statistic.
"""
function Base.summary(io::IO, A::UnivariateStatistic{T, K, I}) where {T, K, I}
    return print(io, "UnivariateStatistic{$T,$K,$I}  with $K moments")
end
#= 
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
end =#

"""
    summary(A::IndependentStatistic)

Return a compact one-line summary of the statistic.
"""
function Base.summary(io::IO, A::IndependentStatistic{T, N, K, I}) where {T, N, K, I}
    return print(io, "$(Base.dims2string(size(A))) IndependentStatistic{$T, $N, $K, $I} with $K moments")
end

# Compact display for REPL and inline use
Base.show(io::IO, A::IndependentStatistic) = Base.summary(io, A)
