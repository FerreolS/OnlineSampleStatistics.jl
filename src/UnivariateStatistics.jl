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
A = UnivariateStatistic(2, [1.0, 2.0, 3.0])
mean(A)  # Calculate the mean (2.0)
var(A)  # Calculate the variance (1.0)
fit!(A, 4.0)  # Add a new data point
```
"""

mutable struct UnivariateStatistic{T, K, I, R} <: OnlineStatsBase.OnlineStat{T}
    weights::I
    rawmoments::R
    @doc "Inner constructor rawmoments can be given as a vector or a NTuple"
    function UnivariateStatistic{T, K, I}(weights::I, rawmoments::Vector{T}) where {T, K, I}
        K == length(rawmoments) || throw(ArgumentError("The length of rawmoments $(length(rawmoments)) must be equal to $K"))
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        !(T <: Complex) || K < 3 || throw(ArgumentError("UnivariateStatistic : $K > 2 not implemented for complex numbers"))

        return new{T, K, I, typeof(rawmoments)}(weights, rawmoments)
    end
    function UnivariateStatistic{T, K, I}(weights::I, rawmoments::NTuple{K, T}) where {T, K, I}
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        !(T <: Complex) || K < 3 || throw(ArgumentError("UnivariateStatistic : $K > 2 not implemented for complex numbers"))
        r = vcat(rawmoments...)
        return new{T, K, I, typeof(r)}(weights, r)
    end

end


"""
    UnivariateStatistic(T::Type, K::Int)

Construct an empty `UnivariateStatistic` with element type `T`, moment order `K`,
and default weight type `Int`.

---

    UnivariateStatistic(T::Type, K::Int, I::Type)

Construct an empty `UnivariateStatistic` with element type `T`, moment order `K`,
and weight type `I`.

    UnivariateStatistic(K::Int, x)

Construct from one sample `x` with moment order `K` and default unit weight.

---

    UnivariateStatistic(K::Int, x, w)

Construct from one sample `x` with one weight `w` and moment order `K`.

---

    UnivariateStatistic(T::Type, K::Int, I::Type, x, w)

Construct from one sample `x` and one weight `w` with explicit value type `T`,
moment order `K`, and weight type `I`.

---

    UnivariateStatistic(T::Type, K::Int, x)

Construct from one sample `x` (converted to `T`) with default weight type `Int`
and unit weight.

---

    UnivariateStatistic(T::Type, K::Int, I::Type, x, w)

Construct from one sample `x` and one weight `w` with element type `T`, moment
order `K`, and weight type `I`.

---

    UnivariateStatistic(T::Type, K::Int, x::AbstractArray)

Construct from an array of samples converted to `T`.

---

    UnivariateStatistic(K::Int, x::AbstractArray)

Construct from an array of samples with moment order `K`, inferring value type
from `x` (integer arrays are promoted to `Float64`).

---

    UnivariateStatistic(K::Int, x::AbstractArray, w::AbstractArray)

Construct from arrays of samples and weights with moment order `K`.

---

    UnivariateStatistic(T::Type, K::Int, I::Type, x::AbstractArray, w::AbstractArray)

Construct from arrays of samples and weights.

---

    UnivariateStatistic(K::Int)

Construct an empty `UnivariateStatistic` of type `Float64` and order `K`.
"""
UnivariateStatistic(T::Type, K::Int) = UnivariateStatistic(T, K, Int)
UnivariateStatistic(K::Int) = UnivariateStatistic(Float64, K)

function UnivariateStatistic(T::Type, K::Int, ::Type{I}) where {I <: Number}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    return build_from_rawmoments(zero(I), zeros(T, K))
end

UnivariateStatistic(T::Type, K::Int, x::Number) = UnivariateStatistic(T, K, Int, x, 1)

UnivariateStatistic(K::Int, x::T) where {T <: Number} = UnivariateStatistic(T, K, x)
UnivariateStatistic(K::Int, x::Integer) = UnivariateStatistic(Float64, K, Float64(x))

UnivariateStatistic(T::Type, K::Int, x::AbstractArray) = UnivariateStatistic(T, K, Int, x)
UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T <: Number} = UnivariateStatistic(T, K, x)
UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T <: Integer} = UnivariateStatistic(Float64, K, x)

function UnivariateStatistic(T::Type{<:Number}, K::Int, ::Type{I}, x::AbstractArray) where {I <: Number}
    A = UnivariateStatistic(T, K, I)
    fit!(A, T.(x))
    return A
end


UnivariateStatistic(K::Int, x::T, w::Number) where {T <: Number} = UnivariateStatistic(T, K, typeof(w), x, w)
UnivariateStatistic(K::Int, x::Integer, w::Number) = UnivariateStatistic(Float64, K, typeof(w), Float64(x), w)

UnivariateStatistic(K::Int, x::AbstractArray{T}, w::AbstractArray{TW}) where {T <: Number, TW <: Number} = UnivariateStatistic(T, K, TW, x, w)
UnivariateStatistic(K::Int, x::AbstractArray{T}, w::AbstractArray{TW}) where {T <: Integer, TW <: Number} = UnivariateStatistic(Float64, K, TW, x, w)

function UnivariateStatistic(T::Type, K::Int, ::Type{I}, x::Number, w::Number) where {I <: Number}
    !(T <: Complex) || K < 3 || throw(ArgumentError("UnivariateStatistic : $K > 2 not implemented for complex numbers"))
    return build_from_rawmoments(I(w), vcat(T(x), zeros(T, K - 1)))
end

function UnivariateStatistic(T::Type, K::Int, ::Type{I}, x::AbstractArray, w::AbstractArray) where {I <: Number}
    size(x) == size(w) || throw(ArgumentError("UnivariateStatistic : size(x) != size(w)"))
    ww = I == Bool ? Int.(w) : I.(w)
    A = UnivariateStatistic(T, K, eltype(ww))
    nonnegative(ww) || throw(ArgumentError("weights can't be negative"))
    fit!(A, T.(x), ww)
    return A
end

UnivariateStatistic(T::Type, K::Int, x::AbstractArray, w::AbstractArray) = UnivariateStatistic(T, K, eltype(w), x, w)

# Backward-compatible wrappers
#= 
@deprecate UnivariateStatistic(x::T, K::Int) where {T <: Number} UnivariateStatistic(T <: AbstractFloat ? T : Float64, K, T <: AbstractFloat ? x : Float64(x))
@deprecate UnivariateStatistic(::Type{T}, x, K::Int) where {T} UnivariateStatistic(T, K, x)
@deprecate UnivariateStatistic(x::AbstractArray{T}, K::Int) where {T <: Real} UnivariateStatistic(Float64, K, Float64.(x))
@deprecate UnivariateStatistic(x::AbstractArray{T}, K::Int) where {T <: Complex} UnivariateStatistic(T, K, x)
@deprecate UnivariateStatistic(x::T, weight::Number, K::Int) where {T <: Number} UnivariateStatistic(T <: AbstractFloat ? T : Float64, K, typeof(weight), T <: AbstractFloat ? x : Float64(x), weight)
@deprecate UnivariateStatistic(T::Type, TW::Type, K::Int) UnivariateStatistic(T, K, TW)
@deprecate UnivariateStatistic(x::AbstractArray{T}, w::AbstractArray, K::Int) where {T <: Real} UnivariateStatistic(Float64, K, Float64.(x), w)
@deprecate UnivariateStatistic(x::AbstractArray{T}, w::AbstractArray{TW}, K::Int) where {T <: Complex, TW <: Number} UnivariateStatistic(T, K, TW, x, w)
 =#
"""
    UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T<:Number}

Constructs a UnivariateStatistic object storing the first `K` moments  from the vector of samples `x`.

"""
function build_from_rawmoments(weights::I, rawmoments::Vector{T}) where {T, I}
    K = length(rawmoments)
    return UnivariateStatistic{T, K, I}(weights, rawmoments)
end

build_from_rawmoments(weights, rawmoments...) = build_from_rawmoments(weights, vcat(rawmoments...))
#UnivariateStatistic{T,K,I}(weights::I, rawmoments...) where {T,K,I} = UnivariateStatistic{T,K,I}(weights, vcat(rawmoments...))

"""
    nonnegative(x)

Check if the input `x` is nonnegative (greater than or equal to zero).
"""
nonnegative(x) = x ≥ 0
nonnegative(x::AbstractArray) = mapreduce(nonnegative, &, x)

"""
    zero(::UnivariateStatistic{T}) where {T<:UnivariateStatistic} 
    zero(::Type{UnivariateStatistic{T,K,I,L}}) where {T,K,I,L}

    return an empty UnivariateStatistic of type `T` with `K` moments and weights of type `I`.
"""
Base.zero(::T) where {T <: UnivariateStatistic} = zero(T)
Base.zero(::Type{UnivariateStatistic{T, K, I, L}}) where {T, K, I, L} = build_from_rawmoments(zero(I), zeros(T, K))

"""
    empty!(A::UnivariateStatistic)

    Reset the UnivariateStatistic `A` to its initial state. 
"""
function Base.empty!(A::UnivariateStatistic)
    A.weights = zero(A.weights)
    A.rawmoments .= zero(eltype(A))
    return A
end


"""
    eltype(::UnivariateStatistic{T}) where {T}

    return the type of the elements in the UnivariateStatistic.
"""

Base.eltype(::UnivariateStatistic{T}) where {T} = T


Base.:(==)(A::UnivariateStatistic{T, K, I}, B::UnivariateStatistic{T, K, I}) where {T, K, I} = A.rawmoments == B.rawmoments && A.weights == B.weights
Base.isapprox(A::UnivariateStatistic{T, K, I}, B::UnivariateStatistic{T, K, I}; kwds...) where {T, K, I} = isapprox(A.rawmoments, B.rawmoments; kwds...) && isapprox(A.weights, B.weights; kwds...)
""" 
    copy(A::UnivariateStatistic)
   Copy  (deepcopy) the UnivariateStatistic `A` to a new object.
"""
Base.copy(A::UnivariateStatistic) = deepcopy(A)


"""
    nobs(A::UnivariateStatistic) 
Return the number of samples  or the sum of weights in a `A`.
"""
StatsBase.nobs(A::UnivariateStatistic) = A.weights
""" 
    weights(A::UnivariateStatistic)
Return the sum of weights in a `A`.
"""
StatsBase.weights(A::UnivariateStatistic{T, K, W}) where {T, K, W <: Number} = A.weights

""""
    order(A::UnivariateStatistic)
Return the number of moments in a `A`.
"""
order(::UnivariateStatistic{T, K}) where {T, K} = K


function get_rawmoments(A::UnivariateStatistic{T, K, I}, k::Int) where {T, K, I}
    k ≤ K || throw(ArgumentError("$k moments are not available for type $(typeof(A))"))
    return A.rawmoments[k]
end
get_rawmoments(A::UnivariateStatistic{T, K, I, R}) where {T, K, I, R} = [get_rawmoments(A, k) for k in 1:K]

"""
    get_moments(A::UnivariateStatistic, k) -> Number

Compute the k-th moment of a UnivariateStatistic `A`. 

"""
get_moments(A::UnivariateStatistic{T, K, I, R}, k) where {T, K, I, R} = ifelse((N = weights(A)) == 0, T(0), get_rawmoments(A, k) / ifelse(k == 1, T(1), T(N)))
get_moments(A::UnivariateStatistic{T, K, I, R}) where {T, K, I, R} = [get_moments(A, k) for k in 1:K]

"""
    mean(A::UnivariateStatistic) 

Compute the sample mean of a  `A` 
"""
Statistics.mean(A::UnivariateStatistic{T, K, I, R}) where {T, K, I, R} = get_moments(A, 1)

"""
    var(A::UnivariateStatistic; corrected=true)

Compute the sample variance of a `A`.  If `corrected` is true, the variance is corrected for bias. 
The unbias variance estimator is only available for an integer number of sample.

"""
function Statistics.var(A::UnivariateStatistic{T, K, W}; corrected = true) where {T, K, W}
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

"""
     std(A::UnivariateStatistic; corrected=true)
Compute the sample standard deviation of a `A`, from its variance (corrected by default).

"""
Statistics.std(A::UnivariateStatistic; corrected = true) = sqrt(var(A; corrected = corrected))


"""
    skewness(A::UnivariateStatistic)
Compute the sample skewness of a `A`. The skewness is defined as the third standardized moment.
"""

function StatsBase.skewness(A::UnivariateStatistic{T, K}) where {T, K}
    3 ≤ K || throw(ArgumentError("third moment is not available for type $(typeof(A))"))
    N = nobs(A)
    N < 2 && return T(NaN)
    cm2, cm3 = A.rawmoments[2:3]
    cm2 == 0 && return T(NaN)
    return cm3 / sqrt(cm2 * cm2 * cm2 / N)
end

"""
    kurtosis(A::UnivariateStatistic)
Compute the sample kurtosis of a `A`. The kurtosis is defined as the fourth standardized moment.
"""

function StatsBase.kurtosis(A::UnivariateStatistic{T, K}) where {T, K}
    4 ≤ K || throw(ArgumentError("fourth moment is not available for type $(typeof(A))"))
    N = nobs(A)
    N < 2 && return T(NaN)
    cm2 = A.rawmoments[2]
    cm2 == 0 && return T(NaN)
    cm4 = A.rawmoments[4]
    return (cm4 / (cm2 * cm2 / N)) - 3.0
end


"""
    fit!(A::UnivariateStatistic{T}, y::T2) where {T, T2}

Pushes a new samples `y` into the UnivariateStatistic `A`.

# Throws
- `ArgumentError`: If the type of elements in `y` is not compatible with the type `T` of `A`.

"""
function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T2}) where {T, T2 <: Number}
    T == eltype(y) || promote_type(T, eltype(y)) == T || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    fit!(A, T.(y))
    return A
end

function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T}) where {T <: Number}
    foreach(x -> _fit!(A, x), y)
    return A
end

function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T2}, w) where {T, T2 <: Number}
    T == eltype(y) || promote_type(T, eltype(y)) == T || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    return fit!(A, T.(y), w)
end

function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T}, w::AbstractArray) where {T <: Number}
    size(y) == size(w) || throw(ArgumentError("UnivariateStatistic : size(y) != size(w)"))
    return fit!(A, zip(y, w))
end

function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T}, w::Real) where {T <: Number}
    if (w == 1)
        foreach(x -> _fit!(A, x), y)
    else
        foreach(x -> _fit!(A, x, w), y)
    end
    return A
end


function fit!(A::UnivariateStatistic{T}, b::T2, w::Real) where {T <: Number, T2 <: Number}
    promote_type(T, T2) == T || throw(ArgumentError("The input type $T2 is not promotable to $T"))
    return _fit!(A, T(b), w)
end

fit!(o::UnivariateStatistic{T}, y::Tuple{T, <:Real}) where {T} = _fit!(o, y...)

@inline increment_weights!(A::UnivariateStatistic, x) = A.weights += x

#= NOT WEIGHTED DATA =#

function _fit!(A::UnivariateStatistic{T, 1}, b::T) where {T <: Number}
    μA, = A.rawmoments
    A.rawmoments[1] += inv(A.weights += 1) * (b - μA)
    return A
end


function _fit!(A::UnivariateStatistic{T, 2}, b::T) where {T <: Number}
    NA = weights(A)
    iN = inv(increment_weights!(A, 1))
    μA, = A.rawmoments

    δBA = (b - μA)
    A.rawmoments[1] += iN * δBA
    A.rawmoments[2] += iN * NA * abs2(δBA)
    return A
end


@generated function _fit!(A::UnivariateStatistic{T, P, Int}, b::T) where {P, T <: Number}
    code = Expr(:block)
    push!(
        code.args, quote
            NA = weights(A)
            N = NA + 1
            iN = inv(increment_weights!(A, 1))
            μA, = A.rawmoments
            δBA = (b - μA)
        end
    )
    for p in P:-1:3
        push!(code.args, :(A.rawmoments[$p] += (NA * (-iN)^$p + (NA * iN)^$p) * δBA^$p))
        for k in 1:(p - 2)
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (A.rawmoments[$p - $k] * (-δBA * iN)^$k)))
        end
    end
    push!(
        code.args, quote
            A.rawmoments[1] += iN * δBA
            A.rawmoments[2] += iN * NA * δBA^2
        end
    )
    push!(code.args, :(return A))
    return code
end


#= WEIGHTED DATA =#

function _fit!(A::UnivariateStatistic{T, 1}, b::T, w::Real) where {T <: Number}
    w == 0 && return A
    A.rawmoments[1] += w * inv(increment_weights!(A, w)) * (b - A.rawmoments[1])
    return A
end


function _fit!(A::UnivariateStatistic{T, 2}, b::T, wb::Real) where {T <: Number}
    wb == 0 && return A
    wa = weights(A)
    μA, _ = A.rawmoments
    iN = inv(increment_weights!(A, wb))
    δBA = (b - μA)
    BoN = -wb * iN * δBA
    AoN = wa * iN * δBA

    A.rawmoments[1] -= BoN
    A.rawmoments[2] += wa * abs2(BoN) + wb * abs2(AoN)

    return A
end

@generated function _fit!(A::UnivariateStatistic{T, P}, b::T, wb::Real) where {P, T <: Number}
    code = Expr(:block)
    push!(
        code.args, quote
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
        end
    )
    for p in P:-1:3
        push!(code.args, :(A.rawmoments[$p] += wa * BoN^$p + wb * AoN^$p))
        for k in 1:(p - 2)
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (A.rawmoments[$p - $k] * BoN^$k)))
        end
    end
    push!(
        code.args, quote
            A.rawmoments[2] += wa * abs2(BoN) + wb * abs2(AoN)
        end
    )
    push!(code.args, :(return A))
    return code
end

"""
    merge(A::UnivariateStatistic, B::UnivariateStatistic)

Merges the statistics from `B` into `A` to a new object. The type of the new object will be promoted.

## Example
```jldoctest
julia> A = UnivariateStatistic(2, [1.0, 0.5])
UnivariateStatistic: n=2 | value=[1.0, 0.25]

julia> B = UnivariateStatistic(2, [2.0, 1.5])
UnivariateStatistic: n=2 | value=[2.0, 0.75]

julia> merge!(A, B)
UnivariateStatistic: n=4 | value=[1.5, 0.75]
```
"""

function Base.merge(A::UnivariateStatistic{T1, 1}, B::UnivariateStatistic{T2, K}) where {T1, T2, K}
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

Merges (inplace) the statistics from `B` into `A` in-place. 

## Example
```julia
A = UnivariateStatistic(2, [1.0, 0.5])
B = UnivariateStatistic(2, [2.0, 1.5])
merge!(A, B) 
A ≈ UnivariateStatistic(2, [1.0, 0.5, 2.0, 1.5])
```
"""
function Base.merge!(A::UnivariateStatistic{T1, 1, I}, B::UnivariateStatistic{T2, K, I}) where {T1, T2, K, I}
    promote_type(T1, T2) == T1 || throw(ArgumentError("The input for $(typeof(A)) is $T1. Found $(eltype(B))."))
    A.weights += B.weights
    A.rawmoments[1] += inv(A.weights) * B.weights * (B.rawmoments[1] - A.rawmoments[1])
    return A
end

function Base.merge!(A::UnivariateStatistic{T1, 2, I}, B::UnivariateStatistic{T2, 2, I}) where {T1, T2, I}
    promote_type(T1, T2) == T1 || throw(ArgumentError("The input for $(typeof(A)) is $T1. Found $(eltype(B))."))
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


@generated function Base.merge!(A::UnivariateStatistic{T, P}, B::UnivariateStatistic{T, M}) where {T, M, P}
    P ≤ M || throw(ArgumentError("The number of moment $M of the second Arguments is less than the first $P."))
    code = Expr(:block)
    push!(
        code.args, quote
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
        end
    )
    for p in P:-1:3
        for k in 1:(p - 2)
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (BoN^$k * A.rawmoments[$p - $k] + AoN^$k * B.rawmoments[$p - $k])))
        end
        push!(code.args, :(A.rawmoments[$p] += B.rawmoments[$p] + wa * BoN^$p + wb * AoN^$p))
    end
    push!(
        code.args, quote
            A.rawmoments[1] -= BoN
            A.rawmoments[2] += MB + wa * abs2(BoN) + wb * abs2(AoN)
        end
    )
    push!(code.args, :(return A))
    return code
end


value(A::UnivariateStatistic) = get_moments(A)
