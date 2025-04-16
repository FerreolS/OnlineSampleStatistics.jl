using ZippedArrays, StructuredArrays, StaticArrays
#= getters on AbstractArrays that are not  IndependentStatistic =#
"""
    get_rawmoments(x::AbstractArray{<:UnivariateStatistic})

Retrieve the raw moments from an array of `UnivariateStatistic` objects. 
Returns an array where each element corresponds to the raw moments of the respective `UnivariateStatistic` in `x`.

"""
function get_rawmoments(A::AbstractArray{<:UnivariateStatistic{T,K},N}) where {T,K,N}
    NTuple{K,Array{T,N}}([get_rawmoments(x, j) for x ∈ A] for j ∈ 1:K)
end

"""

    weights(x::AbstractArray{UnivariateStatistic})

Retrieves the weights from an array of `UnivariateStatistic` objects.
"""
weights(x::AbstractArray{<:UnivariateStatistic}) = map(x -> x.weights, x)

"""
    get_rawmoments(x::AbstractArray{UnivariateStatistic}, k::Int)

Retrieve the `k`-th raw moments from an array of `UnivariateStatistic` objects. 
Returns an array where each element corresponds to the `k`-th raw moment of the respective `UnivariateStatistic` in `x`.

"""
get_rawmoments(x::AbstractArray{<:UnivariateStatistic}, k::Int) = map(y -> get_rawmoments(y, k), x)
"""

    StatsBase.nobs(x::AbstractArray{UnivariateStatistic{T,K,Int}}) where {T,K}

Retrieve the number of observations (`nobs`) from an array of `UnivariateStatistic` objects. 
Returns an array where each element corresponds to the number of observations of the respective `UnivariateStatistic` in `x`.

"""

StatsBase.nobs(x::AbstractArray{<:UnivariateStatistic}) = map(x -> nobs(x), x)

"""
    get_moments(x::AbstractArray{UnivariateStatistic}, k::Int)

Retrieve the `k`-th moments from an array of `UnivariateStatistic` objects. 
Returns an array where each element corresponds to the `k`-th moment of the respective `UnivariateStatistic` in `x`.

"""
get_moments(x::AbstractArray{<:UnivariateStatistic}, k::Int) = map(y -> get_moments(y, k), x)

order(::AbstractArray{<:UnivariateStatistic{T,K}}) where {T,K} = K

#===  IndependentStatistic ===#

IndependentStatistic{T,N,K,W} = ZippedArray{UnivariateStatistic{T,K,W},N,K2,I,A} where {I,A,K2}


#= constructor for IndependentStatistic =#
IndependentStatistic(sz::NTuple{N,Int}, K::Int) where {N} = IndependentStatistic(Float64, sz, K)

function IndependentStatistic(::Type{T}, sz::NTuple{N,Int}, K::Int) where {T,N}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    rawmoments = (zeros(T, sz) for _ in 1:K)
    ZippedArray{UnivariateStatistic{T,K,Int}}(MutableUniformArray(0, sz...), rawmoments...)
end

function IndependentStatistic(::Type{T}, sz::NTuple{N,Int}, ::Type{TW}, K::Int) where {TW,T,N}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    rawmoments = (zeros(T, sz) for _ in 1:K)
    if TW == Bool
        weights = zeros(Int, sz)
    else
        weights = zeros(TW, sz)
    end
    ZippedArray{UnivariateStatistic{T,K,eltype(weights)}}(weights, rawmoments...)
end

function IndependentStatistic(x::AbstractArray{T,N}, w::AbstractArray{TW,N}, K::Int; dims=nothing) where {TW,T,N}
    size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))
    if dims === nothing
        A = IndependentStatistic(T, size(x), TW, K)
        _push!(A, x, w)
    else
        maximum(dims) ≤ N || throw(ArgumentError("IndependentStatistic : $(maximum(dims)) > $N"))
        sz = vcat(size(x)...)
        sz[vcat(dims...)] .= 1
        A = IndependentStatistic(T, NTuple{N,Int}(sz), TW, K)
        #foreach((y, z) -> push!(A, reshape(y, sz...), reshape(z, sz...)), zip(eachslice(x; dims=dims), eachslice(w; dims=dims)))
        for (y, z) ∈ zip(eachslice(x; dims=dims), eachslice(w; dims=dims))
            _push!(A, reshape(y, sz...), reshape(z, sz...))
        end
    end
    return A
end



function IndependentStatistic(x::AbstractArray{T,N}, K::Int; dims=nothing) where {T,N}
    if dims === nothing
        A = IndependentStatistic(T, size(x), K)
        _push!(A, x)
    else
        maximum(dims) ≤ N || throw(ArgumentError("IndependentStatistic : $(maximum(dims)) > $N"))
        sz = vcat(size(x)...)
        sz[vcat(dims...)] .= 1
        A = IndependentStatistic(T, NTuple{N,Int}(sz), K)
        foreach(y -> _push!(A, reshape(y, sz...), 1), eachslice(x; dims=dims))
    end
    return A
end


#= getters on IndependentStatistic =#

get_rawmoments(x::IndependentStatistic{T,N,K,I}) where {T,N,K,I} = @inbounds x.args[2:K+1]
weights(x::IndependentStatistic) = @inbounds x.args[1]
get_rawmoments(x::IndependentStatistic{T,N,K,I}, k::Int) where {T,N,K,I} = @inbounds x.args[2:K+1][k]
order(::IndependentStatistic{T,N,K}) where {T,N,K} = K

#= statistic functions =#
StatsBase.nobs(x::IndependentStatistic) = @inbounds x.args[1]

Statistics.mean(A::IndependentStatistic) = get_rawmoments(A, 1)

function Statistics.var(A::IndependentStatistic{T,N,K,I}; corrected=true) where {T,N,K,I}
    2 ≤ K || throw(ArgumentError("second moment is not available for type $(typeof(A))"))
    W = nobs(A)
    if corrected
        return get_rawmoments(A, 2) ./ (W .- 1)
    else
        return get_rawmoments(A, 2) ./ W
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





function Base.push!(A::IndependentStatistic{T,N}, x::AbstractArray{T,N2}) where {T,N,N2}
    N2 ≥ N || throw(ArgumentError("push! : N2 < N"))
    dims = NTuple{N2 - N,Int}((ndims(A)+1):ndims(x))
    for y ∈ eachslice(x; dims=dims)
        _push!(A, reshape(y, size(A)))
    end
    return A
end

function Base.push!(A::IndependentStatistic{T,N}, x::AbstractArray{T,N}) where {T,N}
    size(A) == size(x) && return _push!(A, x)

    szx = SVector{N,Int}(size(x)...)
    szA = SVector{N,Int}(size(A)...)

    sgltidx = findall(szA .!= 1)
    singleton = ones(MVector{N,Bool})
    singleton[sgltidx] .= false

    szA[sgltidx] == szx[sgltidx] || throw(DimensionMismatch("push! : size(A) incompatible with size(x)"))
    slc = NTuple{sum(singleton),Int}((1:N)[singleton])
    for y ∈ eachslice(x; dims=slc)
        _push!(A, reshape(y, szA...))
    end
    return A
end

@inline function increment!(A::MutableUniformArray, x::Real)
    StructuredArrays.setvalue!(A, StructuredArrays.value(A) + x)
    return A
end
@inline increment!(A::AbstractArray, x) = A .+= x
@inline increment_weights!(A::IndependentStatistic, x) = increment!(weights(A), x)
#= NOT WEIGHTED DATA =#

function _push!(A::IndependentStatistic{T,D,1}, b::AbstractArray{T,D}) where {T<:Number,D}
    wa = weights(A)
    m1 = get_rawmoments(A, 1)
    if wa isa MutableUniformArray
        iN = inv.(increment_weights!(A, T(1)))
        @inbounds @simd for i in eachindex(m1, b)
            δBA = (b[i] - m1[i])
            m1[i] += iN[i] * δBA
        end
    else
        @inbounds @simd for i in eachindex(m1, b)
            wa[i] += T(1)
            δBA = (b[i] - m1[i])
            m1[i] += inv(wa[i]) * δBA
        end
    end

    return A
end



function _push!(A::IndependentStatistic{T,D,2}, b::AbstractArray{T,D}) where {T<:Number,D}
    wa = weights(A)
    m1 = get_rawmoments(A, 1)
    m2 = get_rawmoments(A, 2)
    if wa isa MutableUniformArray
        wa = copy(wa)
        iN = inv.(increment_weights!(A, T(1)))
        @inbounds @simd for i in eachindex(m1, m2, b)
            δBA = (b[i] - m1[i])
            m1[i] += iN[i] * δBA
            m2[i] += iN[i] * wa[i] * abs2(δBA)
        end
    else
        @inbounds @simd for i in eachindex(m1, m2, b)
            δBA = (b[i] - m1[i])
            wa[i] += 1
            iN = inv(wa[i])
            m1[i] += iN * δBA
            m2[i] += iN * wa[i] * abs2(δBA)

        end
    end
    return A
end


@generated function _push!(A::IndependentStatistic{T,D,P}, b::AbstractArray{T,D}) where {D,P,T<:Number}
    code = Expr(:block)
    push!(code.args, quote
        wa = copy(weights(A))
        iN = inv.(increment_weights!(A, 1))
        δBA = @. (b - $get_rawmoments(A, 1))
        @. $get_rawmoments(A, 1) += iN * δBA
        if P == 1
            return A
        end
    end)
    for p in P:-1:3
        push!(code.args, :(@. $get_rawmoments(A, $p) += (wa * (-iN)^$p + (wa * iN)^$p) * δBA^$p))
        for k in 1:(p-2)
            push!(code.args, :(@. $get_rawmoments(A, $p) += binomial($p, $k) * (@. $get_rawmoments(A, $p - $k) * (-δBA * iN)^$k)))
        end
    end
    push!(code.args, quote
        @. $get_rawmoments(A, 2) += iN * wa * abs2(δBA)
    end)
    push!(code.args, :(return A))
    return code
end

#=  WEIGHTED DATA =#
Base.push!(A::IndependentStatistic{T}, x::AbstractArray{T2}, w) where {T,T2} = push!(A, T.(x), w)

function Base.push!(A::IndependentStatistic{T,N}, x::AbstractArray{T,N2}, w::W) where {T,N,N2,W<:Union{Real,AbstractArray{<:Real,N2}}}
    N2 ≥ N || throw(ArgumentError("push! : N2 < N"))
    if W <: AbstractArray
        size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))
    end

    dims = NTuple{N2 - N,Int}((ndims(A)+1):ndims(x))
    for (y, z) ∈ zip(eachslice(x; dims=dims), eachslice(w; dims=dims))
        _push!(A, reshape(y, size(A)), reshape(z, size(A)))
    end
    return A
end

function Base.push!(A::IndependentStatistic{T,N,K}, x::AbstractArray{T,N}, w::W) where {T,N,K,W<:Union{Real,AbstractArray{<:Real,N}}}
    if W <: AbstractArray
        size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))
    end
    szx = SVector{N,Int}(size(x)...)
    szA = SVector{N,Int}(size(A)...)

    sgltidx = findall(szA .!= 1)
    singleton = ones(MVector{N,Bool})
    singleton[sgltidx] .= false

    szA[sgltidx] == szx[sgltidx] || throw(DimensionMismatch("push! : size(A) incompatible with size(x)"))
    slc = NTuple{sum(singleton),Int}((1:N)[singleton])
    if W <: AbstractArray
        for (y, z) ∈ zip(eachslice(x; dims=slc), eachslice(w; dims=slc))
            _push!(A, reshape(y, szA...), reshape(z, szA...))
        end
    else

        for y ∈ eachslice(x; dims=slc)
            _push!(A, reshape(y, szA...), w)
        end
    end
    return A
end

function _push!(A::IndependentStatistic{T,D,1}, b::AbstractArray{T,D}, wb::W) where {D,T<:Real,W<:Union{Real,AbstractArray{<:Real,D}}}

    wa = weights(A)
    m1 = get_rawmoments(A, 1)

    if wa isa MutableUniformArray
        W <: Real || throw(ArgumentError("MutableUniformArray not supported for non scalar weights $W"))
        iN = inv.(increment_weights!(A, wb))

        @inbounds @simd for i in eachindex(m1, b)
            δBA = (b[i] - m1[i])
            m1[i] += iN[i] * δBA * wb
        end
    else
        if W <: AbstractArray
            size(b) == size(wb) || throw(ArgumentError("IndependentStatistic : size(b) != size(wb)"))
            @inbounds @simd for i in eachindex(m1, b, wb)
                wa[i] += wb[i]
                δBA = (b[i] - m1[i])
                m1[i] += inv(wa[i]) * δBA * wb[i]
            end
        else
            @inbounds @simd for i in eachindex(m1, b)
                wa[i] += wb
                δBA = (b[i] - m1[i])
                m1[i] += inv(wa[i]) * δBA * wb
            end
        end
    end
    return A
end

@generated function _push!(A::I, b::AbstractArray{T,D}, wb::W) where {P,D,T<:Real,I<:IndependentStatistic{T,D,P},W<:Union{Real,AbstractArray{<:Real,D}}}
    code = Expr(:block)
    if W <: AbstractArray
        push!(code.args, quote
            size(b) == size(wb) || throw(ArgumentError("IndependentStatistic : size(b) != size(wb)"))
        end)
        WBI = :(wb[i])
    else
        push!(code.args, quote
            wb < 0 && throw(ArgumentError("weights can't be negative"))
            wb == 0 && return A
        end)
        WBI = :(wb)
    end
    push!(code.args, quote
        P > 1 || throw(ArgumentError("P must be greater than 1"))
        nonnegative(wb) || throw(ArgumentError("weights can't be negative"))
        wb == 0 && return A
        wa = weights(A)
        m1 = get_rawmoments(A, 1)
    end)
    if I.parameters[5].parameters[1] <: MutableUniformArray
        W == Real && push!(code.args, :(throw(ArgumentError("MutableUniformArray not supported for non scalar weights"))))
        push!(code.args, :(MWAI = first(wa)))
        push!(code.args, :(MNI = inv.(first(increment_weights!(A, wb)))))
        NI = :(MNI)
        WAI = :(wai = MWAI)
        waupdate = :()
    else
        WAI = :(wai = wa[i])
        NI = :(ifelse($WBI == 0, 0, inv(wa[i] + $WBI)))
        waupdate = :(wa[i] += $WBI)
    end
    push!(code.args, quote
        @inbounds @simd for i in eachindex(m1, b)
            iN = $NI
            $WAI
            wbi = $WBI
            δBA = (b[i] - m1[i])
            BoN = -$WBI * iN * δBA
            m1[i] -= BoN
            AoN = wai * iN * δBA
            $(MomentsExpression(P))
            get_rawmoments(A, 2)[i] += wai * abs2(BoN) + wbi * abs2(AoN)
            $waupdate
        end
        return A
    end)

    return code
end

function MomentsExpression(P::Int)
    code = Expr(:block)
    if P ≤ 2
        return code
    end
    for p in P:-1:3
        push!(code.args, :(get_rawmoments(A, $p)[i] += wai * BoN^$p + wbi * AoN^$p))
        for k in 1:(p-2)
            push!(code.args, :(get_rawmoments(A, $p)[i] += binomial($p, $k) * get_rawmoments(A, $p - $k)[i] * BoN^$k))
        end
    end
    return code
end
