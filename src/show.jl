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
    println(io, "UnivariateStatistic{$T, $K, $I}")
    println(io, "  n: $n")
    
    if n == 0
        println(io, "  (empty)")
        return
    end
    
    # Always show mean
    μ = mean(A)
    print(io, "  μ: ", sprint(show, μ))
    
    # Show variance and std if K ≥ 2
    if K ≥ 2
        σ² = var(A; corrected = false)
        σ = std(A)
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
    
    println(io)
end

# Compact display for REPL and inline use
function Base.show(io::IO, A::UnivariateStatistic{T, K, I}) where {T, K, I}
    n = nobs(A)
    print(io, "UnivariateStatistic{$T,$K}(n=$n)")
end

# Summary method for detailed statistics
"""
    summary(A::UnivariateStatistic)

Print a detailed summary of the statistics including all available computed moments.
"""
function Base.summary(io::IO, A::UnivariateStatistic{T, K, I}) where {T, K, I}
    println(io, "UnivariateStatistic Summary")
    println(io, "==========================")
    
    n = nobs(A)
    println(io, "Observations: $n")
    println(io, "Element type: $T")
    println(io, "Weight type:  $I")
    println(io, "Moments:      $K")
    println(io)
    
    if n == 0
        println(io, "No data.")
        return
    end
    
    println(io, "Statistics:")
    println(io, "  mean          = $(mean(A))")
    
    if K ≥ 2
        println(io, "  variance      = $(var(A; corrected=false))")
        println(io, "  std           = $(std(A))")
    end
    
    if K ≥ 3
        println(io, "  skewness      = $(skewness(A))")
    end
    
    if K ≥ 4
        println(io, "  kurtosis      = $(kurtosis(A))")
    end
    
    println(io)
    println(io, "Raw moments: $(A.rawmoments)")
end

# Pretty-print for IndependentStatistic
function Base.show(io::IO, ::MIME"text/plain", A::IndependentStatistic{T, N, K, I}) where {T, N, K, I}
    sz = size(A)
    wts = weights(A)
    is_uniform = wts isa MutableUniformArray
    
    if is_uniform
        println(io, "IndependentStatistic{$T, $N, $K}")
        println(io, "  size: $sz")
        println(io, "  moments: $K")
        println(io, "  nobs: $(StructuredArrays.value(wts))")
    else
        println(io, "Weighted IndependentStatistic{$T, $N, $K}")
        println(io, "  size: $sz")
        println(io, "  moments: $K")
        
        # Compute weight statistics with reduced precision
        wts_flat = vec(wts)
        wts_min = round(minimum(wts_flat); sigdigits=3)
        wts_max = round(maximum(wts_flat); sigdigits=3)
        wts_median = round(Statistics.median(wts_flat); sigdigits=3)
        println(io, "  weights: min=$wts_min, median=$wts_median, max=$wts_max")
    end
    
    println(io)
end

function Base.show(io::IO, A::IndependentStatistic{T, N, K, I}) where {T, N, K, I}
    sz = size(A)
    print(io, "IndependentStatistic{$T,$N,$K}(size=$sz)")
end
