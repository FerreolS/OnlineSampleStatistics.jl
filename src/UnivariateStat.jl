using Statistics, StatsBase, OnlineStatsBase
"""
    UnivariateStatistic{T,K,I}

A mutable struct for managing univariate statistics, such as mean, variance, skewness, and kurtosis.
This struct supports online (incremental) updates to the statistics.

## Fields
- `weights::I`: the total weight or count (when I==Int) of the data.
- `rawmoments::Vector{T}`: A vector of size K containing the raw moments of the data.

## Example
```julia
A = UnivariateStatistic([1.0, 2.0, 3.0],2,)
mean(A)  # Calculate the mean (2.0)
var(A)  # Calculate the variance (1.0)
push!(A, 4.0)  # Add a new data point
```
"""

mutable struct UnivariateStatistic{T,K,I} <: OnlineStatsBase.OnlineStat{T}
    weights::I
    rawmoments::Vector{T}
    function UnivariateStatistic{T,K,I}(weights::I, rawmoments::Vector{T}) where {T,K,I}
        K == length(rawmoments) || throw(ArgumentError("The length of rawmoments $(length(rawmoments)) must be equal to $K"))
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        new{T,K,I}(weights, rawmoments)
    end
    function UnivariateStatistic{T,K,I}(weights::I, rawmoments::NTuple{K,T}) where {T,K,I}
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        new{T,K,I}(weights, vcat(rawmoments...))
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
    UnivariateStatistic(x::T,K::Int) where {T<:Number}

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

    UnivariateStatistic(K::Int)

Constructs an empty `UnivariateStatistic` object of type Float64 with `K` moments.

---

    UnivariateStatistic(::Type{T}, x, K::Int) where {T}

Constructs a `UnivariateStatistic` object  of type `T` with a single sample `x`.

# Arguments
- `T::Type`: The type to which the elements of `x` will be converted.
- `K::Int`: The number of moments to store.
- `x`: The sample to be converted to type `T` and passed to the constructor.
"""
function UnivariateStatistic(weights::I, rawmoments::Vector{T}) where {T,I}
    K = length(rawmoments)
    UnivariateStatistic{T,K,I}(weights, rawmoments)
end

UnivariateStatistic{T,K,I}(weights::I, rawmoments...) where {T,K,I} = UnivariateStatistic{T,K,I}(weights, vcat(rawmoments...))

UnivariateStatistic(x::T, K::Int) where {T<:AbstractFloat} = UnivariateStatistic(x, 1, K)
UnivariateStatistic(x::T, K::Int) where {T<:Number} = UnivariateStatistic(Float64(x), 1, K)
function UnivariateStatistic(x::T, weight::Number, K::Int) where {T<:Number}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    UnivariateStatistic(weight, vcat(x, zeros(T, K - 1)))
end
UnivariateStatistic(T::Type, K::Int) = UnivariateStatistic(T, Int, K)
UnivariateStatistic(T::Type, TW::Type, K::Int) = UnivariateStatistic(zero(TW), zeros(T, K))
UnivariateStatistic(K::Int) = UnivariateStatistic(Float64, K)
UnivariateStatistic(::Type{T}, x, K::Int) where {T} = UnivariateStatistic(T.(x), K)
"""
    UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T<:Number}

Constructs a UnivariateStatistic object storing the first `K` moments  from the vector of samples `x`.

# Arguments
- `K::Int`: The number of moments 
- `x::AbstractArray{T}`: array of samples where `T` is a subtype of `Number`.
"""
UnivariateStatistic(x::AbstractArray{T}, K::Int) where {T<:Number} = UnivariateStatistic(Float64.(x), 1, K)
function UnivariateStatistic(x::AbstractArray{T}, K::Int) where {T<:Union{AbstractFloat,Complex}}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    !(T <: Complex) || K < 3 || throw(ArgumentError("UnivariateStatistic : $K > 2 not implemented for complex numbers"))
    A = UnivariateStatistic(T, K)
    push!(A, x)
    return A
end

UnivariateStatistic(x::AbstractArray{T}, w::AbstractArray, K::Int) where {T<:Number} = UnivariateStatistic(Float64.(x), w, K)

function UnivariateStatistic(x::AbstractArray{T}, w::AbstractArray{TW}, K::Int) where {T<:Union{AbstractFloat,Complex},TW<:Number}

    size(x) == size(w) || throw(ArgumentError("UnivariateStatistic : size(x) != size(w)"))
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    !(T <: Complex) || K < 3 || throw(ArgumentError("UnivariateStatistic : $K > 2 not implemented for complex numbers"))
    if eltype(w) == Bool
        w = Int.(w)
        A = UnivariateStatistic(T, Int, K)
    else
        A = UnivariateStatistic(T, TW, K)
    end
    nonnegative(w) || throw(ArgumentError("weights can't be negative"))
    push!(A, x, w)
    return A
end

"""
    nonnegative(x)

Check if the input `x` is nonnegative (greater than or equal to zero).
"""
nonnegative(x) = x ≥ 0
nonnegative(x::AbstractArray) = all(x .>= 0)

Base.zero(::T) where {T<:UnivariateStatistic} = zero(T)
Base.zero(::Type{UnivariateStatistic{T,K,I}}) where {T,K,I} = UnivariateStatistic(zero(I), zeros(T, K))
Base.eltype(::UnivariateStatistic{T}) where {T} = T

Base.:(==)(A::UnivariateStatistic{T,K,I}, B::UnivariateStatistic{T,K,I}) where {T,K,I} = A.rawmoments == B.rawmoments && A.weights == B.weights
Base.isapprox(A::UnivariateStatistic{T,K,I}, B::UnivariateStatistic{T,K,I}; kwds...) where {T,K,I} = isapprox(A.rawmoments, B.rawmoments; kwds...) && isapprox(A.weights, B.weights; kwds...)
Base.copy(A::UnivariateStatistic) = deepcopy(A)



StatsBase.nobs(A::UnivariateStatistic) = A.weights
weights(A::UnivariateStatistic{T,K,W}) where {T,K,W<:Number} = A.weights


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
get_moments(A::UnivariateStatistic{T,K}) where {T,K} = [get_moments(A, k) for k in 1:K]

"""
    Statistics.mean(A::UnivariateStatistic{T,K,I}) where {T,K,I}

Compute the mean of a UnivariateStatistic `A` 
# See also
- `get_moments`
"""
Statistics.mean(A::UnivariateStatistic{T,K,I}) where {T,K,I} = get_moments(A, 1)

function Statistics.var(A::UnivariateStatistic{T,K,W}; corrected=true) where {T,K,W}
    N = nobs(A)
    N == 0 && return T(NaN)
    if corrected
        if W <: Integer
            return get_rawmoments(A, 2) / (N - 1)
        else
            @warn("The number of samples is not an integer. The variance is not corrected.")
            return get_rawmoments(A, 2) / N
        end
    else
        return get_rawmoments(A, 2) / N
    end
end


function StatsBase.skewness(A::UnivariateStatistic{T,K,Int}) where {T,K}
    3 ≤ K || throw(ArgumentError("third moment is not available for type $(typeof(A))"))
    N = nobs(A)
    N < 2 && return T(NaN)
    cm2, cm3 = A.rawmoments[2:3]
    cm2 == 0 && return T(NaN)
    return cm3 / sqrt(cm2 * cm2 * cm2 / N)
end

function StatsBase.kurtosis(A::UnivariateStatistic{T,K,Int}) where {T,K}
    3 ≤ K || throw(ArgumentError("third moment is not available for type $(typeof(A))"))
    N = nobs(A)
    N < 2 && return T(NaN)
    cm2 = A.rawmoments[2]
    cm2 == 0 && return T(NaN)
    cm4 = A.rawmoments[4]
    return (cm4 / (cm2 * cm2 / N)) - 3.0
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
function Base.push!(A::UnivariateStatistic{T}, y::AbstractArray{T2}) where {T,T2}
    T == eltype(y) || promote_type(T, eltype(y)) == T || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    for x ∈ y
        A = push!(A, x)
    end
    return A
end

Base.push!(A::UnivariateStatistic, b::Number) = push!(A, b, 1)

function Base.push!(A::UnivariateStatistic{T}, y::AbstractArray{T2}, w::AbstractArray) where {T,T2}
    T == eltype(y) || promote_type(T, eltype(y)) == T || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    for (x, v) ∈ zip(y, w)
        A = push!(A, x, v)
    end
    return A
end

function Base.push!(A::UnivariateStatistic{T}, b::T2, w::Real) where {T,T2<:Number}
    promote_type(T, T2) == T || throw(ArgumentError("The input type $T2 is not promotable to $T"))
    push!(A, T(b), w)
end



@inline increment_weights!(A::UnivariateStatistic, x) = A.weights += x

@generated function Base.push!(A::UnivariateStatistic{T,P}, b::T, wb::Real) where {P,T<:Number}
    code = Expr(:block)
    push!(code.args, quote
        wb == 0 && return A
        wa = weights(A)
        iN = inv(increment_weights!(A, wb))
        δBA = (b - A.rawmoments[1])
        BoN = -wb * iN * δBA
        A.rawmoments[1] -= BoN
        if P == 1
            return A
        end
        AoN = wa * iN * δBA
    end)
    for p in P:-1:3
        push!(code.args, :(A.rawmoments[$p] += wa * BoN^$p + wb * AoN^$p))
        for k in 1:(p-2)
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (A.rawmoments[$p-$k] * BoN^$k)))
        end
    end
    push!(code.args, quote
        A.rawmoments[2] += wa * abs2(BoN) + wb * abs2(AoN)
    end)
    push!(code.args, :(return A))
    return code
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
function Base.merge(A::UnivariateStatistic{T1,1}, B::UnivariateStatistic{T2,K}) where {T1,T2,K}
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
    (wb = weights(B)) == 0 && return A
    wa = weights(A)
    iN = inv(increment_weights!(A, wb))
    μA = A.rawmoments[1]
    μB, MB = B.rawmoments[1:2]
    δBA = (μB - μA)
    BoN = -wb * iN * δBA
    AoN = wa * iN * δBA

    # A.rawmoments[2] += MB + inv(N) * wa * wb * δBA^2
    A.rawmoments[1] -= BoN
    A.rawmoments[2] += MB + wa * abs2(BoN) + wb * abs2(AoN)
    return A
end


@generated function Base.merge!(A::UnivariateStatistic{T,P}, B::UnivariateStatistic{T,M}) where {T,M,P}
    P ≤ M || throw(ArgumentError("The number of moment $M of the second Arguments is less than the first $P."))
    code = Expr(:block)
    push!(code.args, quote
        (wb = weights(B)) == 0 && return A
        wa = weights(A)
        iN = inv(increment_weights!(A, wb))
        μA = A.rawmoments[1]
        μB, MB = B.rawmoments[1:2]
        δBA = (μB - μA)
        BoN = -wb * iN * δBA
        AoN = wa * iN * δBA

        BoN = -iN * wb * δBA
        AoN = iN * wa * δBA
    end)
    for p in P:-1:3
        for k in 1:(p-2)
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (BoN^$k * A.rawmoments[$p-$k] + AoN^$k * B.rawmoments[$p-$k])))
        end
        push!(code.args, :(A.rawmoments[$p] += B.rawmoments[$p] + wa * BoN^$p + wb * AoN^$p))
    end
    push!(code.args, quote
        A.rawmoments[1] -= BoN
        A.rawmoments[2] += MB + wa * abs2(BoN) + wb * abs2(AoN)
    end)
    push!(code.args, :(return A))
    return code
end


#Base.merge!(A::UnivariateStatistic, x::Number) = push!(A, x)
OnlineStatsBase._fit!(A::UnivariateStatistic, x::Number) = push!(A, x)

OnlineStatsBase._fit!(A::UnivariateStatistic, x::Number, w::Real) = push!(A, x, w)

#OnlineStatsBase._fit!(A::UnivariateStatistic, x::Base.Iterators.Zip) = push!(A, x.is[1], x.is[2])
#OnlineStatsBase._fit!(A::UnivariateStatistic, x::AbstractArray) = push!(A, x)

OnlineStatsBase.value(A::UnivariateStatistic) = get_moments(A)

function Base.empty!(A::UnivariateStatistic)
    A.weights = zero(A.weights)
    A.rawmoments .= zero(eltype(A))
    return A
end


OnlineStatsBase._merge!(A::UnivariateStatistic, B::UnivariateStatistic) = merge!(A, B)

#= Overloading fit! for UnivariateStatistic to make Transducers working for weighted data =#

function OnlineStatsBase.fit!(o::UnivariateStatistic{I}, y::Iterators.Zip{<:Tuple{AbstractArray{T},<:AbstractArray{<:Real}}}) where {I,T}
    I == T || error("The input for $(name(o,false,false)) is $I. Found $T.")
    for (yi, wi) in y
        OnlineStatsBase.fit!(o, yi, wi)
    end
    o
end
OnlineStatsBase.fit!(o::UnivariateStatistic{T}, y::T, w::Real) where {T} = (OnlineStatsBase._fit!(o, y, w); return o)
OnlineStatsBase.fit!(o::UnivariateStatistic{T}, y::Tuple{T,<:Real}) where {T} = OnlineStatsBase._fit!(o, y...)