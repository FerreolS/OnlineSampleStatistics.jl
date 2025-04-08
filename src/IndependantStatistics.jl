using ZippedArrays, StructuredArrays

#= getters on AbstractArrays that are not  IndependentStatistic =#
get_rawmoments(x::AbstractArray{UnivariateStatistic}) = map(x -> x.rawmoments, x)
get_weights(x::AbstractArray{UnivariateStatistic}) = map(x -> x.weights, x)
weights(x::AbstractArray{UnivariateStatistic}) = map(x -> x.weights, x)
get_rawmoments(x::AbstractArray{UnivariateStatistic}, k::Int) = map(y -> get_rawmoments(y, k), x)
StatsBase.nobs(x::AbstractArray{UnivariateStatistic{T,K,Int}}) where {T,K} = map(x -> nobs(x), x)
get_moments(x::AbstractArray{UnivariateStatistic}, k::Int) = map(y -> get_moments(y, k), x)

#===  IndependentStatistic ===#

IndependentStatistic{T,N,K,W} = ZippedArray{UnivariateStatistic{T,K,W},N,K2,I,A} where {I,A,K2}


#= constructor for IndependentStatistic =#
IndependentStatistic(K::Int, sz::NTuple{N,Int}) where {N} = IndependentStatistic(Float64, K, sz)

function IndependentStatistic(::Type{T}, K::Int, sz::NTuple{N,Int}) where {T,N}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    rawmoments = (zeros(T, sz) for _ in 1:K)
    ZippedArray{UnivariateStatistic{T,K,Int}}(MutableUniformArray(0, sz...), rawmoments...)
end

function IndependentStatistic(::Type{T}, K::Int, sz::NTuple{N,Int}, ::Type{TW}) where {TW,T,N}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    rawmoments = (zeros(T, sz) for _ in 1:K)
    if TW == Bool
        weights = zeros(Int, sz)
    else
        weights = zeros(TW, sz)
    end
    ZippedArray{UnivariateStatistic{T,K,eltype(weights)}}(weights, rawmoments...)
end

function IndependentStatistic(K::Int, x::AbstractArray{T,N}, w::AbstractArray{TW,N}; dims=nothing) where {TW,T,N}
    size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))
    if dims === nothing
        A = IndependentStatistic(T, K, size(x), TW)
        push!(A, x, w)
    else
        maximum(dims) ≤ N || throw(ArgumentError("IndependentStatistic : $(maximum(dims)) > $N"))
        sz = vcat(size(x)...)
        sz[vcat(dims...)] .= 1
        A = IndependentStatistic(T, K, NTuple{N,Int}(sz), TW)
        #foreach((y, z) -> push!(A, reshape(y, sz...), reshape(z, sz...)), zip(eachslice(x; dims=dims), eachslice(w; dims=dims)))
        for (y, z) ∈ zip(eachslice(x; dims=dims), eachslice(w; dims=dims))
            _push!(A, reshape(y, sz...), reshape(z, sz...))
        end
    end
    return A
end



function IndependentStatistic(K::Int, x::AbstractArray{T,N}; dims=nothing) where {T,N}
    if dims === nothing
        A = IndependentStatistic(T, K, size(x))
        push!(A, x, 1)
    else
        maximum(dims) ≥ N || throw(ArgumentError("IndependentStatistic : $(maximum(dims)) > $N"))
        sz = vcat(size(x)...)
        sz[vcat(dims...)] .= 1
        A = IndependentStatistic(T, K, NTuple{N,Int}(sz))
        foreach(y -> _push!(A, reshape(y, sz...), 1), eachslice(x; dims=dims))
    end
    return A
end


#= getters on IndependentStatistic =#

get_rawmoments(x::IndependentStatistic{T,N,K,I}) where {T,N,K,I} = @inbounds x.args[2:K+1]
get_weights(x::IndependentStatistic) = @inbounds x.args[1]
weights(x::IndependentStatistic) = @inbounds x.args[1]
get_rawmoments(x::IndependentStatistic, k::Int) = @inbounds x.args[1+k]

#= statistic functions =#
StatsBase.nobs(x::IndependentStatistic) = @inbounds x.args[1]

Statistics.mean(A::IndependentStatistic) = get_rawmoments(A, 1)

function Statistics.var(A::IndependentStatistic{T,N,K}; corrected=true) where {T,N,K}
    2 ≤ K || throw(ArgumentError("second moment is not available for type $(typeof(A))"))
    W = nobs(A)
    if corrected
        return @. get_rawmoments(A, 2) / (W - 1)
    else
        return @. get_rawmoments(A, 2) / W
    end
end


function StatsBase.skewness(A::IndependentStatistic{T,N,K}) where {T,N,K}
    3 ≤ K || throw(ArgumentError("third moment is not available for type $(typeof(A))"))
    W = nobs(A)
    cm2 = get_rawmoments(A, 2)
    cm3 = get_rawmoments(A, 3)
    return @. cm3 / sqrt(cm2 * cm2 * cm2 / W)
end

function StatsBase.kurtosis(A::IndependentStatistic{T,N,K}) where {T,N,K}
    4 ≤ K || throw(ArgumentError("fourth moment is not available for type $(typeof(A))"))
    W = nobs(A)
    cm2 = get_rawmoments(A, 2)
    cm4 = get_rawmoments(A, 4)
    return @. (cm4 / (cm2 * cm2 / W)) - 3.0
end




Base.push!(A::IndependentStatistic, x) = push!(A, x, 1)
Base.push!(A::IndependentStatistic{T}, x::AbstractArray{T2}, w) where {T,T2} = push!(A, T.(x), w)

function Base.push!(A::IndependentStatistic{T,N,K,W}, x::AbstractArray{T,N2}, w::AbstractArray{<:Real,N2}) where {T,N,N2,K,W}
    N2 ≥ N || throw(ArgumentError("push! : N2 < N :  $N2 < $N; $(size(x)) and $(size(A))"))
    size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))

    szx = vcat(size(x)...)
    szA = vcat(size(A)...)

    sgltidx = findall(szA .!= 1)
    singleton = trues(N2)
    singleton[sgltidx] .= false

    szA[sgltidx] == szx[sgltidx] || throw(ArgumentError("push! : size(A)=$(size(A)) != $(size(x))"))
    slc = NTuple{sum(singleton),Int}((1:N2)[singleton])

    for (y, z) ∈ zip(eachslice(x; dims=slc, drop=(length(sgltidx) == length(szA))), eachslice(w; dims=slc, drop=(length(sgltidx) == length(szA))))
        _push!(A, reshape(y, szA...), reshape(z, szA...))
    end
    return A
end

function Base.push!(A::IndependentStatistic{T,N}, x::AbstractArray{T,N2}, w::Real) where {T,N,N2}
    N2 ≥ N || throw(ArgumentError("push! : N2 < N"))

    szx = vcat(size(x)...)
    szA = vcat(size(A)...)

    sgltidx = findall(szA .!= 1)
    singleton = trues(N2)
    singleton[sgltidx] .= false

    szA[sgltidx] == szx[sgltidx] || throw(ArgumentError("push! : size(A)=$(size(A)) != $(size(x))"))
    slc = NTuple{sum(singleton),Int}((1:N2)[singleton])
    for y ∈ eachslice(x; dims=slc, drop=(length(sgltidx) == length(szA)))
        _push!(A, reshape(y, szA...), w)
    end
    return A
end



@inline function increment!(A::MutableUniformArray, x::Real)
    StructuredArrays.setvalue!(A, StructuredArrays.value(A) + x)
    return A
end
@inline increment!(A::AbstractArray, x) = A .+= x
@inline increment_weights!(A::IndependentStatistic, x) = increment!(weights(A), x)

@generated function _push!(A::IndependentStatistic{T,D,P}, b::AbstractArray{T,D}, wb::W) where {P,D,T<:Number,W<:Union{Real,AbstractArray{<:Number,D}}}
    code = Expr(:block)
    push!(code.args, quote
        wa = copy(weights(A))
        iN = inv.(increment_weights!(A, wb))
        @. iN[!isfinite(iN)] = zero(T)
        δBA = @. (b - $get_rawmoments(A, 1))
        BoN = @. -wb * iN * δBA
        @. $get_rawmoments(A, 1) -= BoN
        if P == 1
            return A
        end
        AoN = @. wa * iN * δBA
    end)
    for p in P:-1:3
        push!(code.args, :(@. $get_rawmoments(A, $p) += wa * BoN^$p + wb * AoN^$p))
        for k in 1:(p-2)
            push!(code.args, :(@. $get_rawmoments(A, $p) += binomial($p, $k) * (@. $get_rawmoments(A, $p - $k) * BoN^$k)))
        end
    end
    push!(code.args, quote
        @. $get_rawmoments(A, 2) += wa * abs2(BoN) + wb * abs2(AoN)
    end)
    push!(code.args, :(return A))
    return code
end
