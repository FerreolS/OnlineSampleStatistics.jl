struct MyStructT{N}
    data::NTuple{N,Int}
end

struct MyStructNT{T}
    data::@NamedTuple{T}
end

struct MyStructA
    data::Vector{Int}
end


export UnivariateStatistic2, update



struct UnivariateStatistic2{T,K,I,R}
    weights::I
    #rawmoments::Vector{T}
    rawmoments::R
    function UnivariateStatistic2{T,K,I}(weights::I, rawmoments) where {T,K,I}
        #  nonnegative(weights) || throw(ArgumentError("weights can't be negative"))
        new{T,K,I,typeof(rawmoments)}(weights, rawmoments)
    end
end

weights(A::UnivariateStatistic2{T,K,W}) where {T,K,W<:Number} = A.weights


function update(A::UnivariateStatistic2{T,1}, b::T) where {T<:Number}
    μA, = A.rawmoments
    NA = weights(A)
    A.rawmoments = A.rawmoments .+ (inv(NA + 1) * (b - μA))
    return UnivariateStatistic2(NA + 1, (μA,))
end


function update(A::UnivariateStatistic2{T,2,I}, b::T) where {T<:Number,I}
    NA = A.weights
    iN = inv(NA + 1)

    μA, Ma = A.rawmoments

    δBA = (b - μA)
    return UnivariateStatistic2{T,2,I}(weights(A) + 1, (μA + iN * δBA, Ma + iN * weights(A) * abs2(δBA)))
end

f(A, x) =
    foreach(x) do b
        A = update(A, b)
    end