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
    rawmoments = (zeros(T, sz) for _ in 1:K)
    ZippedArray{UnivariateStatistic{T,K,Int}}(MutableUniformArray(0, sz...), rawmoments...)
end

function IndependentStatistic(K::Int, x::AbstractArray{T,N}; dims=nothing) where {T,N}
    if dims === nothing
        A = IndependentStatistic(T, K, size(x))
        push!(A, x)
    else
        sz = vcat(size(x)...)
        sz[vcat(dims...)] .= 1
        A = IndependentStatistic(T, K, NTuple{N,Int}(sz))
        foreach(y -> push!(A, reshape(y, sz...)), eachslice(x; dims=dims))
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





Base.push!(A::IndependentStatistic{T}, x::AbstractArray{T2}) where {T,T2} = push!(A, T.(x))

function Base.push!(A::IndependentStatistic{T,N,K,W}, x::AbstractArray{T,N2}) where {T,N,N2,K,W}
    N2 ≥ N || throw(ArgumentError("push! : N2 < N"))
    szx = vcat(size(x)...)
    szA = vcat(size(A)...)

    sgltidx = findall(szA .!= 1)
    singleton = trues(N2)
    singleton[sgltidx] .= false

    szA[sgltidx] == szx[sgltidx] || throw(ArgumentError("push! : size(A)=$(size(A)) != $(size(x))"))
    # mapslices(y -> _push!(A, reshape(y, szA...)), x; dims=(1:ndims(x))[.!singleton])
    slc = NTuple{sum(singleton),Int}((1:N2)[singleton])
    #foreach(y -> _push!(A, reshape(y, szA...)), eachslice(x; dims=slc))
    foreach(y -> push!(A, y), eachslice(x; dims=slc, drop=(length(sgltidx) == length(szA))))
    return A
end



function add!(A::MutableUniformArray, x::Number)
    StructuredArrays.setvalue!(A, StructuredArrays.value(A) + x)
    return A
end
add!(A::AbstractArray, x) = A .+= x


function _push!(A::IndependentStatistic{T,N,1}, b::AbstractArray{T,N}) where {T,N}
    w = get_weights(A)
    add!(w, 1)
    @. $get_rawmoments(A, 1) += inv(w) * (b - $get_rawmoments(A, 1))
    return A
end

function _push!(A::IndependentStatistic{T,D,2}, b::AbstractArray{T,D}) where {T,D}
    N = add!(weights(A), 1)
    NA = N .- 1
    δBA = (b .- get_rawmoments(A, 1))
    iN = inv.(N)
    @. $get_rawmoments(A, 1) += iN * δBA
    @. $get_rawmoments(A, 2) += iN * NA * δBA^2
    return A
end

@generated function Base.push!(A::IndependentStatistic{T,D,P}, b::AbstractArray{T,D}) where {P,D,T<:Number}
    code = Expr(:block)
    if P < 3
        push!(code.args, :(_push!(A, b)))
        return code
    end
    push!(code.args, quote
        N = add!(weights(A), 1)
        NA = N .- 1
        δBA = (b .- get_rawmoments(A, 1))
        iN = inv.(N)
    end)
    for p in P:-1:3
        push!(code.args, :(@. $get_rawmoments(A, $p) += (NA * (-iN)^$p + (NA * iN)^$p) * δBA^$p))
        for k in 1:(p-2)
            push!(code.args, :(@. $get_rawmoments(A, $p) += binomial($p, $k) * (@. $get_rawmoments(A, $p - $k) * (-δBA * iN)^$k)))
        end
    end
    push!(code.args, quote
        @. $get_rawmoments(A, 1) += iN * δBA
        @. $get_rawmoments(A, 2) += iN * NA * δBA^2
    end)
    push!(code.args, :(return A))
    return code
end
