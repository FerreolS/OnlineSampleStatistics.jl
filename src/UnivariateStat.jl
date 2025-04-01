using StatsAPI, Statistics, StatsBase

mutable struct UnivariateStatistic{T,K,I}
    rawmoments::Vector{T}
    weights::I
    function UnivariateStatistic(rawmoments::Vector{T}, weights::I) where {T,I}
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        K = length(rawmoments)
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
UnivariateStatistic(K::Int, x::T) where {T<:Number} = UnivariateStatistic(vcat(x, zeros(T, K - 1)), 1)
UnivariateStatistic(K::Int, T::Type) = UnivariateStatistic(zeros(T, K), 0)
UnivariateStatistic(::Type{T}, K::Int, x) where {T} = UnivariateStatistic(K, T.(x))

function UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T<:Number}
    A = UnivariateStatistic(K, T)
    push!(A, x)
    return A
end

nonnegative(x) = x ≥ 0

Base.zero(::T) where {T<:UnivariateStatistic} = zero(T)
Base.zero(::Type{UnivariateStatistic{T,K,I}}) where {T,K,I} = UnivariateStatistic(zeros(T, K), zero(I))
Base.eltype(::UnivariateStatistic{T}) where {T} = T

Base.:(==)(A::UnivariateStatistic{T,K,I}, B::UnivariateStatistic{T,K,I}) where {T,K,I} = A.rawmoments == B.rawmoments && A.weights == B.weights
Base.copy(A::UnivariateStatistic) = deepcopy(A)

StatsAPI.nobs(A::UnivariateStatistic{T,K,Int}) where {T,K} = A.weights
StatsAPI.weights(A::UnivariateStatistic) = A.weights


function get_rawmoments(A::UnivariateStatistic{T,K,I}, k::Int) where {T,K,I}
    k ≤ K || throw(ArgumentError("$k moments are not available for type $(typeof(A))"))
    return A.rawmoments[k]
end

get_moments(A::UnivariateStatistic{T,K,I}, k) where {T,K,I} = ifelse((N = weights(A)) == 0, 0, get_rawmoments(A, k) / ifelse(k == 1, 1, N))

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

function Base.merge(A::UnivariateStatistic{T1,1,I}, B::UnivariateStatistic{T2,K,I}) where {T1,T2,K,I}
    T = promote_type(T1, T2)
    if T == T1
        C = copy(A)
        merge!(C, B)
        return C
    elseif T == T2
        C = copy(B)
        merge!(C, A)
        return C
    else
        error("The input for $(typeof(A)) is $T. Found $T2.")
    end
end


function Base.merge!(A::UnivariateStatistic{T1,1,I}, B::UnivariateStatistic{T2,K,I}) where {T1,T2,K,I}
    promote_type(T1, T2) == T1 || error("The input for $(typeof(A)) is $T. Found $(eltype(B)).")
    A.weights += B.weights
    A.rawmoments[1] += inv(A.weights) * B.weights * (B.rawmoments[1] - A.rawmoments[1])
    return A
end


function Base.merge!(A::UnivariateStatistic{T1,2,I}, B::UnivariateStatistic{T2,2,I}) where {T1,T2,I}
    promote_type(T1, T2) == T1 || error("The input for $(typeof(A)) is $T. Found $(eltype(B)).")
    NA = weights(A)
    (NB = weights(B)) == 0 && return A
    N = NA + NB
    A.weights = N
    μA, MA = A.rawmoments
    μB, MB = B.rawmoments
    δAB = (μB - μA)
    A.rawmoments[1] += inv(N) * NB * δAB
    A.rawmoments[2] += MB + inv(N) * NA * NB * δAB^2
    return A
end

function Base.push!(A::UnivariateStatistic{T}, y::T2) where {T,T2}
    T == eltype(y) || promote_type(T, eltype(y)) == T || error("The input for $(typeof(A)) is $T. Found $T2.")
    for x ∈ y
        A = push!(A, x)
    end
    return A
end

function Base.push!(A::UnivariateStatistic{T}, b::T2) where {T,T2<:Number}
    promote_type(T, T2) == T || error("The input type $T2 is not promotable to $T")
    push!(A, T(b))
end


function Base.push!(A::UnivariateStatistic{T,1}, b::T) where {T<:Number}
    A.rawmoments[1] += inv(A.weights += 1) * (b - A.rawmoments[1])
    return A
end

function Base.push!(A::UnivariateStatistic{T,2}, b::T) where {T<:Number}
    NA = weights(A)
    N = NA + 1
    A.weights = N
    μA, MA = A.rawmoments

    δAB = (b - μA)
    A.rawmoments[1] += inv(N) * δAB

    A.rawmoments[2] += inv(N) * NA * δAB^2
    return A
end
