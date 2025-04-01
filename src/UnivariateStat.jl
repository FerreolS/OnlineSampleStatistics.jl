using StatsAPI, Statistics, StatsBase
"""
    UnivariateStatistic{T,K,I}

A mutable struct for managing univariate statistics, such as mean, variance, skewness, and kurtosis.
This struct supports online (incremental) updates to the statistics.

## Fields
- `rawmoments::Vector{T}`: A vector of size K containing the raw moments of the data.
- `weights::I`: the total weight or count (when I==Int) of the data.

## Example
```julia
A = UnivariateStatistic(2, [1.0, 2.0, 3.0])
mean(A)  # Calculate the mean (2.0)
var(A)  # Calculate the variance (1.0)
push!(A, 4.0)  # Add a new data point
```
"""

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
"""
    UnivariateStatistic(K::Int, x::T) where {T<:Number}

Constructs a `UnivariateStatistic` object of type `T` with `K` moments from a single sample `x`. 
The first moment (the mean) is then `x` and the remaining moments are zeros of type `T`. 
The weight counts the number of sample set to `1`.

# Arguments
- `K::Int`: The number of moments to store.
- `x::T`: The first sample.

---

    UnivariateStatistic(K::Int, T::Type)

Constructs an empty `UnivariateStatistic` object of type `T` with `K` moments.
# Arguments
- `K::Int`: The number of moments to store.
- `T::Type`: The type of the elements in the moments.

---

    UnivariateStatistic(::Type{T}, K::Int, x) where {T}

Constructs a `UnivariateStatistic` object  of type `T` with a single sample `x`.

# Arguments
- `T::Type`: The type to which the elements of `x` will be converted.
- `K::Int`: The number of moments to store.
- `x`: The sample to be converted to type `T` and passed to the constructor.
"""
UnivariateStatistic(K::Int, x::T) where {T<:Number} = UnivariateStatistic(vcat(x, zeros(T, K - 1)), 1)
UnivariateStatistic(K::Int, T::Type) = UnivariateStatistic(zeros(T, K), 0)
UnivariateStatistic(::Type{T}, K::Int, x) where {T} = UnivariateStatistic(K, T.(x))
"""
    UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T<:Number}

Constructs a UnivariateStatistic object storing the first `K` moments  from the vector of samples `x`.

# Arguments
- `K::Int`: The number of moments 
- `x::AbstractArray{T}`: array of samples where `T` is a subtype of `Number`.
"""

function UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T<:Number}
    A = UnivariateStatistic(K, T)
    push!(A, x)
    return A
end

"""
    nonnegative(x)

Check if the input `x` is nonnegative (greater than or equal to zero).
"""
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

"""
    get_moments(A::UnivariateStatistic{T,K,I}, k) -> Number

Compute the k-th moment of a UnivariateStatistic `A`. 

# Arguments
- `A::UnivariateStatistic{T,K,I}`: The UnivariateStatistic object containing the samples statistic.
- `k::Int`: The order of the moment to compute.

# Returns
- The k-th moment of the statistic. If it is empty, the function returns `0`. 
"""
get_moments(A::UnivariateStatistic{T,K,I}, k) where {T,K,I} = ifelse((N = weights(A)) == 0, 0, get_rawmoments(A, k) / ifelse(k == 1, 1, N))

"""
    Statistics.mean(A::UnivariateStatistic{T,K,I}) where {T,K,I}

Compute the mean of a UnivariateStatistic `A` 
# See also
- `get_moments`
"""
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
    return sqrt(get_rawmoments(A, 3) / N)
end

function StatsBase.kurtosis(A::UnivariateStatistic{T,K,Int}) where {T,K}
    N = nobs(A)
    N < 2 && return T(NaN)
    return get_rawmoments(A, 4) / N
end

"""
    merge(A::UnivariateStatistic, B::UnivariateStatistic)

Merges the statistics from `B` into `A` to a new object. The type of the new object will be promoted.

## Arguments
- `A::UnivariateStatistic`: The destination statistic to be updated.
- `B::UnivariateStatistic`: The source statistic to merge into `A`.

## Example
```julia
A = UnivariateStatistic(2, [1.0, 0.5])
B = UnivariateStatistic(2, [2.0, 1.5])
merge!(A, B)
```
"""
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
        throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    end
end

"""
    merge!(A::UnivariateStatistic, B::UnivariateStatistic)

Merges (inplace) the statistics from `B` into `A` in-place. The weights and raw moments of `A` are updated
to reflect the combined statistics of both `A` and `B`.

## Arguments
- `A::UnivariateStatistic`: The destination statistic to be updated.
- `B::UnivariateStatistic`: The source statistic to merge into `A`.

## Example
```julia
A = UnivariateStatistic(2, [1.0, 0.5])
B = UnivariateStatistic(2, [2.0, 1.5])
merge!(A, B)
```
"""
function Base.merge!(A::UnivariateStatistic{T1,1,I}, B::UnivariateStatistic{T2,K,I}) where {T1,T2,K,I}
    promote_type(T1, T2) == T1 || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $(eltype(B))."))
    A.weights += B.weights
    A.rawmoments[1] += inv(A.weights) * B.weights * (B.rawmoments[1] - A.rawmoments[1])
    return A
end

function Base.merge!(A::UnivariateStatistic{T1,2,I}, B::UnivariateStatistic{T2,2,I}) where {T1,T2,I}
    promote_type(T1, T2) == T1 || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $(eltype(B))."))
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

"""
    Base.push!(A::UnivariateStatistic{T}, y::T2) where {T, T2}

Pushes a new samples `y` into the UnivariateStatistic `A`.

# Arguments
- `A::UnivariateStatistic{T}`: The UnivariateStatistic object to which elements will be added. The type parameter `T` specifies the expected element type.
- `y::T2`: the sample or array of sample pushed to `A`. The type `T2` must either match `T` or be promotable to `T`.

# Throws
- `ArgumentError`: If the type of elements in `y` is not compatible with the type `T` of `A`.

"""
function Base.push!(A::UnivariateStatistic{T}, y::T2) where {T,T2}
    T == eltype(y) || promote_type(T, eltype(y)) == T || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    for x ∈ y
        A = push!(A, x)
    end
    return A
end

function Base.push!(A::UnivariateStatistic{T}, b::T2) where {T,T2<:Number}
    promote_type(T, T2) == T || throw(ArgumentError("The input type $T2 is not promotable to $T"))
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



@generated function Base.merge!(A::UnivariateStatistic{T,N,Int}, B::UnivariateStatistic{T,M,Int}) where {T,M,N}
    N < M || throw(ArgumentError("The number of moment $M of the second Arguments is less than the first $N."))
    code = Expr(:block)
    push!(code.args, quote
        NA = weights(A)
        (NB = weights(B)) == 0 && return A
        N = NA + NB
        A.weights = N
        μA, MA = A.rawmoments[1:2]
        μB, MB = B.rawmoments[1:2]
        δAB = (μB - μA)
        P1 = -inv(N) * NB * δAB
        P2 = inv(N) * NA * δAB
        A.rawmoments[1] -= P1
        A.rawmoments[2] += MB + P2 * NB * δAB
    end)
    for p in 3:N
        for k in 0:p
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (P1^$k * A.rawmoments[$p-$k] + P2^$k * B.rawmoments[$p-$k])))
        end
    end
    push!(code.args, :(return A))
    return code
end


@generated function Base.push!(A::UnivariateStatistic{T,P,Int}, b::T) where {P,T<:Number}
    code = Expr(:block)
    push!(code.args, quote
        NA = weights(A)
        N = NA + 1
        A.weights = N
        μA, MA = A.rawmoments[1:2]
        δAB = (b - μA)
        A.rawmoments[1] += inv(N) * δAB
        A.rawmoments[2] += inv(N) * NA * δAB^2
    end)
    for p in 3:P
        #push!(code.args, :(A.rawmoments[$p] += ((inv(-N))^$p * (N - 1) + (inv(N) * (N - 1))^$p) * δAB))
        push!(code.args, :(A.rawmoments[$p] += ((N - 1) / (-N)^$p + ((N - 1) / N)^$p) * δAB^$p))
        for k in 1:(p-2)
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (A.rawmoments[$p-$k] * (inv(-N))^$k * δAB)^$k))
        end
    end
    push!(code.args, :(return A))
    return code
end
