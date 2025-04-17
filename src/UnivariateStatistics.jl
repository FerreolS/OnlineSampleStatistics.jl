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
fit!(A, 4.0)  # Add a new data point
```
"""

mutable struct UnivariateStatistic{T,K,I,R} <: OnlineStatsBase.OnlineStat{T}
    weights::I
    rawmoments::R
    @doc "Inner constructor rawmoments can be given as a vector or a NTuple"
    function UnivariateStatistic{T,K,I}(weights::I, rawmoments::Vector{T}) where {T,K,I}
        K == length(rawmoments) || throw(ArgumentError("The length of rawmoments $(length(rawmoments)) must be equal to $K"))
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        new{T,K,I,typeof(rawmoments)}(weights, rawmoments)
    end
    function UnivariateStatistic{T,K,I}(weights::I, rawmoments::NTuple{K,T}) where {T,K,I}
        nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        r = vcat(rawmoments...)
        new{T,K,I,typeof(r)}(weights, r)
    end

end


"""
    UnivariateStatistic(x::T,K::Int) where {T<:Number}

Constructs a `UnivariateStatistic` object of type `T` with `K` moments from a single sample `x`. 
The first moment (the mean) is then `x` and the remaining moments are zeros of type `T`. 
The weight counts the number of sample set to `1`.

---

    UnivariateStatistic(K::Int, T::Type)

Constructs an empty `UnivariateStatistic` object of type `T` with `K` moments.

---

    UnivariateStatistic(K::Int)

Constructs an empty `UnivariateStatistic` object of type Float64 with `K` moments.

---

    UnivariateStatistic(::Type{T}, x, K::Int) where {T}

Constructs a `UnivariateStatistic` object  of type `T` with a single sample `x`. 
`x`will be converted to type `T` is needed.

"""
UnivariateStatistic(x::T, K::Int) where {T<:AbstractFloat} = UnivariateStatistic(x, 1, K)
UnivariateStatistic(x::T, K::Int) where {T<:Number} = UnivariateStatistic(Float64(x), 1, K)
UnivariateStatistic(T::Type, K::Int) = UnivariateStatistic(T, Int, K)
UnivariateStatistic(K::Int) = UnivariateStatistic(Float64, K)
UnivariateStatistic(::Type{T}, x, K::Int) where {T} = UnivariateStatistic(T.(x), K)
UnivariateStatistic(x::AbstractArray{T}, K::Int) where {T<:Number} = UnivariateStatistic(Float64.(x), 1, K)

"""
    UnivariateStatistic(K::Int, x::AbstractArray{T}) where {T<:Number}

Constructs a UnivariateStatistic object storing the first `K` moments  from the vector of samples `x`.

"""
function UnivariateStatistic(x::AbstractArray{T}, K::Int) where {T<:Union{AbstractFloat,Complex}}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    !(T <: Complex) || K < 3 || throw(ArgumentError("UnivariateStatistic : $K > 2 not implemented for complex numbers"))
    A = UnivariateStatistic(T, K)
    fit!(A, x)
    return A
end

"""
    UnivariateStatistic(weights::I, rawmoments::Vector{T}) where {T,I}

Constructs a `UnivariateStatistic` object of type `T` with moments given in `rawmoments`.
The `weights` is the sum of weights of the samples that were used to construct `rawmoments`
.
"""

function UnivariateStatistic(weights::I, rawmoments::Vector{T}) where {T,I}
    K = length(rawmoments)
    UnivariateStatistic{T,K,I}(weights, rawmoments)
end

UnivariateStatistic{T,K,I}(weights::I, rawmoments...) where {T,K,I} = UnivariateStatistic{T,K,I}(weights, rawmoments)
#UnivariateStatistic{T,K,I}(weights::I, rawmoments...) where {T,K,I} = UnivariateStatistic{T,K,I}(weights, vcat(rawmoments...))

"""
    UnivariateStatistic(x::T, weight::Number, K::Int) where {T<:Number}

Constructs a `UnivariateStatistic` object with `K` moments of type `T` from a single sample `x` and a weight `weight`.
The first moment (the mean) is then `x` and the remaining moments are zeros of type `T`.
"""
function UnivariateStatistic(x::T, weight::Number, K::Int) where {T<:Number}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    UnivariateStatistic(weight, vcat(x, zeros(T, K - 1)))
end
"""
    UnivariateStatistic(T::Type, TW::Type, K::Int)
Constructs an empty `UnivariateStatistic` object of type `T` with `K` moments and weights of type `TW`.
"""
UnivariateStatistic(T::Type, TW::Type, K::Int) = UnivariateStatistic(zero(TW), zeros(T, K))

"""
    UnivariateStatistic(x::AbstractArray{T}, w::AbstractArray, K::Int) where {T<:Number}      
Constructs a `UnivariateStatistic` object of type `T` with `K` moments from a vector of samples `x` and a vector of weights `w`.
"""
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
    fit!(A, x, w)
    return A
end

"""
    nonnegative(x)

Check if the input `x` is nonnegative (greater than or equal to zero).
"""
nonnegative(x) = x ≥ 0
nonnegative(x::AbstractArray) = all(x .>= 0)

"""
    zero(::UnivariateStatistic{T}) where {T<:UnivariateStatistic} 
    zero(::Type{UnivariateStatistic{T,K,I,L}}) where {T,K,I,L}

    return an empty UnivariateStatistic of type `T` with `K` moments and weights of type `I`.
"""
Base.zero(::T) where {T<:UnivariateStatistic} = zero(T)
Base.zero(::Type{UnivariateStatistic{T,K,I,L}}) where {T,K,I,L} = UnivariateStatistic(zero(I), zeros(T, K))

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


Base.:(==)(A::UnivariateStatistic{T,K,I}, B::UnivariateStatistic{T,K,I}) where {T,K,I} = A.rawmoments == B.rawmoments && A.weights == B.weights
Base.isapprox(A::UnivariateStatistic{T,K,I}, B::UnivariateStatistic{T,K,I}; kwds...) where {T,K,I} = isapprox(A.rawmoments, B.rawmoments; kwds...) && isapprox(A.weights, B.weights; kwds...)
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
StatsBase.weights(A::UnivariateStatistic{T,K,W}) where {T,K,W<:Number} = A.weights

""""
    order(A::UnivariateStatistic)
Return the number of moments in a `A`.
"""
order(::UnivariateStatistic{T,K}) where {T,K} = K


function get_rawmoments(A::UnivariateStatistic{T,K,I}, k::Int) where {T,K,I}
    k ≤ K || throw(ArgumentError("$k moments are not available for type $(typeof(A))"))
    return A.rawmoments[k]
end

"""
    get_moments(A::UnivariateStatistic, k) -> Number

Compute the k-th moment of a UnivariateStatistic `A`. 

"""
get_moments(A::UnivariateStatistic{T,K,I,R}, k) where {T,K,I,R} = ifelse((N = weights(A)) == 0, T(0), get_rawmoments(A, k) / ifelse(k == 1, T(1), T(N)))
get_moments(A::UnivariateStatistic{T,K,I,R}) where {T,K,I,R} = [get_moments(A, k) for k in 1:K]

"""
    mean(A::UnivariateStatistic) 

Compute the sample mean of a  `A` 
"""
Statistics.mean(A::UnivariateStatistic{T,K,I,R}) where {T,K,I,R} = get_moments(A, 1)

"""
    var(A::UnivariateStatistic; corrected=true)

Compute the sample variance of a `A`.  If `corrected` is true, the variance is corrected for bias. 
The unbias variance estimator is only available for an integer number of sample.

"""
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

"""
    std(A::UnivariateStatistic)
Compute the uncorrected sample standard deviation of a `A`. 

"""
Statistics.std(A::UnivariateStatistic) = sqrt(var(A))

"""
    skewness(A::UnivariateStatistic)
Compute the sample skewness of a `A`. The skewness is defined as the third standardized moment.
"""

function StatsBase.skewness(A::UnivariateStatistic{T,K,Int}) where {T,K}
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
    fit!(A::UnivariateStatistic{T}, y::T2) where {T, T2}

Pushes a new samples `y` into the UnivariateStatistic `A`.

# Throws
- `ArgumentError`: If the type of elements in `y` is not compatible with the type `T` of `A`.

"""
function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T2}) where {T,T2}
    T == eltype(y) || promote_type(T, eltype(y)) == T || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    fit!(A, T.(y))
    return A
end

fit!(A::UnivariateStatistic{T}, y::AbstractArray{T}) where {T} = foreach(x -> _fit!(A, x), y)

function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T2}, w) where {T,T2}
    T == eltype(y) || promote_type(T, eltype(y)) == T || throw(ArgumentError("The input for $(typeof(A)) is $T. Found $T2."))
    fit!(A, T.(y), w)
end

function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T}, w::AbstractArray) where {T}
    for (x, v) ∈ zip(y, w)
        _fit!(A, x, v)
    end
    return A
end

function fit!(A::UnivariateStatistic{T}, y::AbstractArray{T}, w::Real) where {T}
    if (w == 1)
        foreach(x -> _fit!(A, x), y)
    else
        foreach(x -> _fit!(A, x, w), y)
    end
    return A
end


function fit!(A::UnivariateStatistic{T}, b::T2, w::Real) where {T<:Number,T2<:Number}
    promote_type(T, T2) == T || throw(ArgumentError("The input type $T2 is not promotable to $T"))
    _fit!(A, T(b), w)
end



@inline increment_weights!(A::UnivariateStatistic, x) = A.weights += x

#= NOT WEIGHTED DATA =#

function _fit!(A::UnivariateStatistic{T,1}, b::T) where {T<:Number}
    μA, = A.rawmoments
    A.rawmoments[1] += inv(A.weights += 1) * (b - μA)
    return A
end


function _fit!(A::UnivariateStatistic{T,2}, b::T) where {T<:Number}
    NA = weights(A)
    iN = inv(increment_weights!(A, 1))
    μA, = A.rawmoments

    δBA = (b - μA)
    A.rawmoments[1] += iN * δBA
    A.rawmoments[2] += iN * NA * abs2(δBA)
    return A
end


@generated function _fit!(A::UnivariateStatistic{T,P,Int}, b::T) where {P,T<:Number}
    code = Expr(:block)
    push!(code.args, quote
        NA = weights(A)
        N = NA + 1
        iN = inv(increment_weights!(A, 1))
        μA, = A.rawmoments
        δBA = (b - μA)
    end)
    for p in P:-1:3
        push!(code.args, :(A.rawmoments[$p] += (NA * (-iN)^$p + (NA * iN)^$p) * δBA^$p))
        for k in 1:(p-2)
            push!(code.args, :(A.rawmoments[$p] += binomial($p, $k) * (A.rawmoments[$p-$k] * (-δBA * iN)^$k)))
        end
    end
    push!(code.args, quote
        A.rawmoments[1] += iN * δBA
        A.rawmoments[2] += iN * NA * δBA^2
    end)
    push!(code.args, :(return A))
    return code
end


#= WEIGHTED DATA =#

function _fit!(A::UnivariateStatistic{T,1}, b::T, w::Real) where {T<:Number}
    w == 0 && return A
    A.rawmoments[1] += w * inv(increment_weights!(A, w)) * (b - A.rawmoments[1])
    return A
end


function _fit!(A::UnivariateStatistic{T,2}, b::T, wb::Real) where {T<:Number}
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

@generated function _fit!(A::UnivariateStatistic{T,P}, b::T, wb::Real) where {P,T<:Number}
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

Merges (inplace) the statistics from `B` into `A` in-place. 

## Example
```julia
A = UnivariateStatistic(2, [1.0, 0.5])
B = UnivariateStatistic(2, [2.0, 1.5])
merge!(A, B) 
A ≈ UnivariateStatistic(2, [1.0, 0.5, 2.0, 1.5])
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


#Base.merge!(A::UnivariateStatistic, x::Number) = fit!(A, x)
#OnlineStatsBase._fit!(A::UnivariateStatistic, x::Number) = fit!(A, x)

#OnlineStatsBase._fit!(A::UnivariateStatistic, x::Number, w::Real) = fit!(A, x, w)

#OnlineStatsBase._fit!(A::UnivariateStatistic, x::Base.Iterators.Zip) = fit!(A, x.is[1], x.is[2])
#OnlineStatsBase._fit!(A::UnivariateStatistic, x::AbstractArray) = fit!(A, x)

value(A::UnivariateStatistic) = get_moments(A)

#OnlineStatsBase._merge!(A::UnivariateStatistic, B::UnivariateStatistic) = merge!(A, B)

#= Overloading fit! for UnivariateStatistic to make Transducers working for weighted data =#

# function OnlineStatsBase.fit!(o::UnivariateStatistic{I}, y::Iterators.Zip{<:Tuple{AbstractArray{T},<:AbstractArray{<:Real}}}) where {I,T}
#     I == T || error("The input for $(name(o,false,false)) is $I. Found $T.")
#     for (yi, wi) in y
#         OnlineStatsBase.fit!(o, yi, wi)
#     end
#     o
# end
#OnlineStatsBase.fit!(o::UnivariateStatistic{T}, y::T, w::Real) where {T} = OnlineStatsBase._fit!(o, y, w)
fit!(o::UnivariateStatistic{T}, y::Tuple{T,<:Real}) where {T} = _fit!(o, y...)