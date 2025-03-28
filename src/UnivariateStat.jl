using StatsAPI, Statistics, StatsBase

struct UnivariateStatistic{T,K,I} # <: Real?
    rawmoments::NTuple{K,T}
    weights::I
    function UnivariateStatistic(rawmoments::NTuple{K,T}, weights::I) where {T,K,I}
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        new{T,K,I}(rawmoments, weights)
    end
end

#= Future traits system 
abstract type WeightTraits end
struct FrequencyWeights <: WeightTraits end
struct AnalyticWeights <: WeightTraits end
struct ProbabilityWeights <: WeightTraits end
WeightTraits(A::UnivariateStatistic) = ProbabilityWeights()
WeightTraits(A::UnivariateStatistic{T,K,Int}) where {T,K} = FrequencyWeights()
 =#
UnivariateStatistic(K::Int, x::T) where {T<:Number} = UnivariateStatistic((x, Tuple(zeros(K - 1))...), 1)
UnivariateStatistic(K::Int, T::Type) = UnivariateStatistic(Tuple(zeros(T, K)), 0)

nonnegative(x) = x ≥ 0

Base.zero(::Type{UnivariateStatistic{T,K,I}}) where {T,K,I} = UnivariateStatistic(Tuple(zeros(T, N)), zero(I))

StatsAPI.nobs(A::UnivariateStatistic{T,K,Int}) where {T,K} = A.weights
StatsAPI.weights(A::UnivariateStatistic) = A.weights


function get_rawmoments(A::UnivariateStatistic{T,K,I}, k::Int) where {T,K,I}
    k ≤ K || throw(ArgumentError("$k moments are not available for type $(typeof(A))"))
    return A.rawmoments[k]
end

get_moments(A::UnivariateStatistic{T,K,I}, k) where {T,K,I} = ifelse((N = nobs(A)) == 0, 0, get_rawmoments(A, k) / N)

Statistics.mean(A::UnivariateStatistic{T,K,I}) where {T,K,I} = get_moments(A, 1)

function Statistics.var(A::UnivariateStatistic{T,K,Int}; corrected=true) where {T,K}
    N = nobs(A)
    N == 0 && return T(NaN)
    if corrected
        return get_rawmoments(A, 2) / (N - 1)
    else
        return get_rawmoments(A, 2) / N
    end
end


function StatsBase.skewness(A::UnivariateStatistic{T,K,Int}) where {T,K}
    N = nobs(A)
    N < 2 && return T(NaN)
    return get_rawmoments(A, 3) / N
end

function StatsBase.kurtosis(A::UnivariateStatistic{T,K,Int}) where {T,K}
    N = nobs(A)
    N < 2 && return T(NaN)
    return get_rawmoments(A, 4) / N
end

function Base.merge(A::UnivariateStatistic{TA,2,Int}, B::UnivariateStatistic{TB,2,Int}) where {TA,TB}


    NA = nobs(A)
    NB = nobs(B)
    N = NA + NB

    μA = get_moments(A, 1)
    μB = get_moments(B, 1)
    δBA = μB - μA
    μ = μA + δBA * NB / N

    MA = get_rawmoments(A, 2)
    MB = get_rawmoments(B, 2)
    M = MA + MB + δBA^2 * NA * NB / N
    return UnivariateStatistic((μ, M), N)
end

function Base.merge(A::UnivariateStatistic{T,2,Int}, b::Number) where {T}
    NA = nobs(A)
    N = NA + 1

    μA = get_moments(A, 1)
    δBA = b - μA
    μ = μA + δBA / N

    MA = get_rawmoments(A, 2)
    M = MA + δBA^2 * NA / N
    return UnivariateStatistic((μ, M), N)
end

function StatsAPI.fit(A::UnivariateStatistic{T,2,Int}, V) where {T}
    for x ∈ V
        A = merge(A, x)
    end
    return A
end