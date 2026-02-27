import ZippedArrays: ZippedArray, ZippedArrays
import StructuredArrays
import StructuredArrays: MutableUniformArray, setvalue!
#= getters on AbstractArrays that are not  IndependentStatistic =#
"""
    get_rawmoments(x::AbstractArray{<:UnivariateStatistic})

Retrieve the raw moments from an array of `UnivariateStatistic` objects. 
Returns an array where each element corresponds to the raw moments of the respective `UnivariateStatistic` in `x`.

"""
function get_rawmoments(A::AbstractArray{<:UnivariateStatistic{T, K}, N}) where {T, K, N}
    return NTuple{K, Array{T, N}}([get_rawmoments(x, j) for x in A] for j in 1:K)
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

order(::AbstractArray{<:UnivariateStatistic{T, K}}) where {T, K} = K

#===  IndependentStatistic ===#

IndependentStatistic{T, N, K, W} = ZippedArray{UnivariateStatistic{T, K, W}, N, K2, I, A} where {I, A, K2}

ZippedArrays.build(::Type{<:UnivariateStatistic}, args...) = build_from_rawmoments(args...)
ZippedArrays.build(::Type{<:UnivariateStatistic}, args::Tuple) = build_from_rawmoments(args...)

#= constructor for IndependentStatistic =#
IndependentStatistic(K::Int, sz::NTuple{N, Int}) where {N} = IndependentStatistic(Float64, K, sz)

function IndependentStatistic(::Type{T}, K::Int, sz::NTuple{N, Int}) where {T, N}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    rawmoments = (zeros(T, sz) for _ in 1:K)
    return ZippedArray{UnivariateStatistic{T, K, Int}}(MutableUniformArray(0, sz...), rawmoments...)
end

function IndependentStatistic(::Type{T}, K::Int, ::Type{TW}, sz::NTuple{N, Int}) where {TW, T, N}
    K > 0 || throw(ArgumentError("Moment of order $K <= 0 undefined"))
    rawmoments = (zeros(T, sz) for _ in 1:K)
    if TW == Bool
        weights = zeros(Int, sz)
    else
        weights = zeros(TW, sz)
    end
    return ZippedArray{UnivariateStatistic{T, K, eltype(weights)}}(weights, rawmoments...)
end

@deprecate IndependentStatistic(sz::NTuple{N, Int}, K::Int) where {N} IndependentStatistic(K, sz)
@deprecate IndependentStatistic(::Type{T}, sz::NTuple{N, Int}, K::Int) where {T, N} IndependentStatistic(T, K, sz)
@deprecate IndependentStatistic(::Type{T}, sz::NTuple{N, Int}, ::Type{TW}, K::Int) where {TW, T, N} IndependentStatistic(T, K, TW, sz)

@inline function _check_fit_compatible_size(szA::NTuple{N, Int}, szx::NTuple{N2, Int}) where {N, N2}
    @inbounds for d in 1:N
        (szA[d] == 1 || szA[d] == szx[d]) ||
            throw(DimensionMismatch("fit! : size(A) incompatible with size(x)"))
    end
    return nothing
end

@generated function _sliced_fit!(::Val{SZ}, A::IndependentStatistic{T, N, K}, x::AbstractArray{T, N2}) where {T, N, N2, K, SZ}
    slc = tuple(tuple(findall(n -> n == 1, SZ)...)..., ((N + 1):N2)...)
    code = Expr(:block)
    push!(
        code.args, :(
            for y in eachslice(x; dims = $slc)
                _fit!(A, reshape(y, SZ))
            end
        )
    )
    push!(code.args, :(return A))
    return code
end

@generated function _sliced_fit!(::Val{SZ}, A::IndependentStatistic{T, N, K}, x::AbstractArray{T, N2}, w::AbstractArray{TW, N2}) where {T, N, N2, TW, K, SZ}
    slc = tuple(tuple(findall(n -> n == 1, SZ)...)..., ((N + 1):N2)...)

    code = Expr(:block)
    push!(
        code.args, :(
            for (y, z) in zip(eachslice(x; dims = $slc), eachslice(w; dims = $slc))
                _fit!(A, reshape(y, SZ), reshape(z, SZ))
            end
        )
    )
    push!(code.args, :(return A))
    return code
end

@generated function _sliced_fit!(::Val{SZ}, A::IndependentStatistic{T, N, K}, x::AbstractArray{T, N2}, w::Number) where {T, N, N2, K, SZ}
    slc = tuple(tuple(findall(n -> n == 1, SZ)...)..., ((N + 1):N2)...)
    code = Expr(:block)
    push!(
        code.args, :(
            for y in eachslice(x; dims = $slc)
                _fit!(A, reshape(y, SZ), T(w))
            end
        )
    )
    push!(code.args, :(return A))
    return code
end

function IndependentStatistic(K::Int, x::AbstractArray{T, N}, w::AbstractArray{TW, N}; dims = nothing) where {TW, T, N}
    size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))
    if dims === nothing
        A = IndependentStatistic(T, K, TW, size(x))
        _fit!(A, x, w)
    else
        maximum(dims) ≤ N || throw(ArgumentError("IndependentStatistic : $(maximum(dims)) > $N"))
        sz = size(x)
        szA = ntuple(i -> ((i ∈ dims) ? 1 : sz[i]), N)
        A = IndependentStatistic(T, K, TW, NTuple{N, Int}(szA))
        _sliced_fit!(Val(szA), A, x, w)

    end
    return A
end


function IndependentStatistic(K::Int, x::AbstractArray{T, N}; dims = nothing) where {T, N}
    if dims === nothing
        A = IndependentStatistic(T, K, size(x))
        _fit!(A, x)
    else
        maximum(dims) ≤ N || throw(ArgumentError("IndependentStatistic : $(maximum(dims)) > $N"))
        sz = size(x)
        szA = ntuple(i -> ((i ∈ dims) ? 1 : sz[i]), N)
        A = IndependentStatistic(T, K, NTuple{N, Int}(szA))
        _sliced_fit!(Val(szA), A, x)
    end
    return A
end

function IndependentStatistic(x::AbstractArray{T, N}, w::AbstractArray{TW, N}, K::Int; dims = nothing) where {TW, T, N}
    Base.depwarn("`IndependentStatistic(x, w, K; dims=...)` is deprecated, use `IndependentStatistic(K, x, w; dims=...)`.", :IndependentStatistic)
    return IndependentStatistic(K, x, w; dims = dims)
end

function IndependentStatistic(x::AbstractArray{T, N}, K::Int; dims = nothing) where {T, N}
    Base.depwarn("`IndependentStatistic(x, K; dims=...)` is deprecated, use `IndependentStatistic(K, x; dims=...)`.", :IndependentStatistic)
    return IndependentStatistic(K, x; dims = dims)
end


#= getters on IndependentStatistic =#

get_rawmoments(x::IndependentStatistic{T, N, K, I}) where {T, N, K, I} = @inbounds x.args[2:(K + 1)]
weights(x::IndependentStatistic) = @inbounds x.args[1]
get_rawmoments(x::IndependentStatistic{T, N, K, I}, k::Int) where {T, N, K, I} = @inbounds x.args[2:(K + 1)][k]
order(::IndependentStatistic{T, N, K}) where {T, N, K} = K

#= statistic functions =#
StatsBase.nobs(x::IndependentStatistic) = @inbounds x.args[1]

Statistics.mean(A::IndependentStatistic) = get_rawmoments(A, 1)

function Statistics.var(A::IndependentStatistic{T, N, K, I}; corrected = true) where {T, N, K, I}
    2 ≤ K || throw(ArgumentError("second moment is not available for type $(typeof(A))"))
    W = nobs(A)
    if corrected
        return get_rawmoments(A, 2) ./ (W .- 1)
    else
        return get_rawmoments(A, 2) ./ W
    end
end


function StatsBase.skewness(A::IndependentStatistic{T, N, K}) where {T, N, K}
    3 ≤ K || throw(ArgumentError("third moment is not available for type $(typeof(A))"))
    W = nobs(A)
    cm2 = get_rawmoments(A, 2)
    cm3 = get_rawmoments(A, 3)
    return @. cm3 / sqrt(cm2 * cm2 * cm2 / W)
end

function StatsBase.kurtosis(A::IndependentStatistic{T, N, K}) where {T, N, K}
    4 ≤ K || throw(ArgumentError("fourth moment is not available for type $(typeof(A))"))
    W = nobs(A)
    cm2 = get_rawmoments(A, 2)
    cm4 = get_rawmoments(A, 4)
    return @. (cm4 / (cm2 * cm2 / W)) - 3.0
end


function fit!(A::IndependentStatistic{T, N}, x::AbstractArray{T, N2}) where {T, N, N2}
    N2 ≥ N || throw(ArgumentError("fit! : N2 < N"))
    _sliced_fit!(Val(size(A)), A, x)
    return A
end

function fit!(A::IndependentStatistic{T, N}, x::AbstractArray{T, N}) where {T, N}
    size(A) == size(x) && return _fit!(A, x)

    szx = size(x)
    szA = size(A)
    _check_fit_compatible_size(szA, szx)
    _sliced_fit!(Val(szA), A, x)
    return A
end

@inline function increment!(A::MutableUniformArray, x::Real)
    setvalue!(A, StructuredArrays.value(A) + x)
    return A
end
@inline increment!(A::AbstractArray, x) = A .+= x
@inline increment_weights!(A::IndependentStatistic, x) = increment!(weights(A), x)
#= NOT WEIGHTED DATA =#

function _fit!(A::IndependentStatistic{T, D, 1}, b::AbstractArray{T, D}) where {T <: Number, D}
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


function _fit!(A::IndependentStatistic{T, D, 2}, b::AbstractArray{T, D}) where {T <: Number, D}
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
            w = wa[i]
            wa[i] += 1
            iN = inv(wa[i])
            m1[i] += iN * δBA
            m2[i] += iN * w * abs2(δBA)
        end
    end
    return A
end


@generated function _fit!(A::IndependentStatistic{T, D, P}, b::AbstractArray{T, D}) where {D, P, T <: Number}
    code = Expr(:block)
    push!(
        code.args, quote
            wa = copy(weights(A))
            iN = inv.(increment_weights!(A, 1))
            δBA = @. (b - $get_rawmoments(A, 1))
            @. $get_rawmoments(A, 1) += iN * δBA
            if P == 1
                return A
            end
        end
    )
    for p in P:-1:3
        push!(code.args, :(@. $get_rawmoments(A, $p) += (wa * (-iN)^$p + (wa * iN)^$p) * δBA^$p))
        for k in 1:(p - 2)
            push!(code.args, :(@. $get_rawmoments(A, $p) += binomial($p, $k) * (@. $get_rawmoments(A, $p - $k) * (-δBA * iN)^$k)))
        end
    end
    push!(
        code.args, quote
            @. $get_rawmoments(A, 2) += iN * wa * abs2(δBA)
        end
    )
    push!(code.args, :(return A))
    return code
end

#=  WEIGHTED DATA =#
fit!(A::IndependentStatistic{T}, x::AbstractArray{T2}, w) where {T, T2} = fit!(A, T.(x), w)

function fit!(A::IndependentStatistic{T, N}, x::AbstractArray{T, N2}, w::W) where {T, N, N2, W <: Union{Real, AbstractArray{<:Real, N2}}}
    N2 ≥ N || throw(ArgumentError("fit! : N2 < N"))
    if W <: AbstractArray
        size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))
    end

    szx = size(x)
    szA = size(A)
    _check_fit_compatible_size(szA, szx)
    _sliced_fit!(Val(szA), A, x, w)
    return A
end

function fit!(A::IndependentStatistic{T, N, K}, x::AbstractArray{T, N}, w::W) where {T, N, K, W <: Union{Real, AbstractArray{<:Real, N}}}
    if W <: AbstractArray
        size(x) == size(w) || throw(ArgumentError("IndependentStatistic : size(x) != size(w)"))
    end

    szx = size(x)
    szA = size(A)
    _check_fit_compatible_size(szA, szx)
    _sliced_fit!(Val(szA), A, x, w)
    return A
end

"""
    fit!(A::IndependentStatistic{T,N,K}, B::IndependentStatistic)

Merge the accumulated statistics in `B` into `A` in-place.
`A` must have compatible element type and size, and its number of moments `K`
must be less than or equal to the one of `B`.
"""
function fit!(A::IndependentStatistic{T, N, K}, B::IndependentStatistic{TB, N, KB}) where {T, TB, N, K, KB}
    K ≤ KB || throw(ArgumentError("fit! : order(B) = $KB is less than order(A) = $K"))
    promote_type(T, TB) == T || throw(ArgumentError("fit! : incompatible element types $TB for $T"))
    size(A) == size(B) || throw(DimensionMismatch("fit! : size(A) != size(B)"))

    wA = weights(A)
    wB = weights(B)
    wA_mutable = wA isa MutableUniformArray
    wB_mutable = wB isa MutableUniformArray
    wA_mutable && !wB_mutable && throw(ArgumentError("fit! : cannot merge non-uniform weights into MutableUniformArray weights"))

    mA = get_rawmoments(A)
    mB = get_rawmoments(B)

    uniform_new_weight = nothing
    @inbounds for i in eachindex(get_rawmoments(A, 1), get_rawmoments(B, 1))
        wa = wA_mutable ? StructuredArrays.value(wA) : wA[i]
        wb = wB_mutable ? StructuredArrays.value(wB) : wB[i]
        wnew = wa + wb

        if wA_mutable
            if uniform_new_weight === nothing
                uniform_new_weight = wnew
            else
                wnew == uniform_new_weight || throw(ArgumentError("fit! : resulting weights are not uniform"))
            end
        else
            wA[i] = wnew
        end

        wb == 0 && continue

        μA = mA[1][i]
        μB = mB[1][i]
        iN = inv(wnew)
        δBA = μB - μA
        BoN = -wb * iN * δBA
        AoN = wa * iN * δBA

        if K ≥ 3
            for p in K:-1:3
                for k in 1:(p - 2)
                    mA[p][i] += binomial(p, k) * (BoN^k * mA[p - k][i] + AoN^k * mB[p - k][i])
                end
                mA[p][i] += mB[p][i] + wa * BoN^p + wb * AoN^p
            end
        end

        mA[1][i] -= BoN
        if K ≥ 2
            mA[2][i] += mB[2][i] + wa * abs2(BoN) + wb * abs2(AoN)
        end
    end

    wA_mutable && setvalue!(wA, uniform_new_weight)
    return A
end

"""
    merge(A::IndependentStatistic, B::IndependentStatistic)

Merge two independent statistics by combining their data.

Creates a copy of statistic `A` and fits it with the data from statistic `B`,
returning a new merged statistic `C`.

# Arguments
- `A::IndependentStatistic`: The first independent statistic
- `B::IndependentStatistic`: The second independent statistic to merge into `A`

# Returns
- `C::IndependentStatistic`: A new statistic that is the result of merging `A` and `B`

# See Also
- [`merge!`](@ref): In-place version of merge

"""
function Base.merge(A::IndependentStatistic, B::IndependentStatistic)
    C = deepcopy(A)
    fit!(C, B)
    return C
end

"""
    merge!(A::IndependentStatistic, B::IndependentStatistic)

Merge the statistics from independent statistic `B` into `A`, updating `A` in-place.

This operation combines the sample statistics from `B` into `A` by fitting `A` with the data 
represented by `B`. After merging, `A` will contain aggregated statistics from both sources.

# Arguments
- `A::IndependentStatistic`: The target statistic object to be updated (modified in-place)
- `B::IndependentStatistic`: The source statistic object to merge into `A`

# Returns
- `A::IndependentStatistic`: The updated statistic object
"""
Base.merge!(A::IndependentStatistic, B::IndependentStatistic) = fit!(A, B)

function _fit!(A::IndependentStatistic{T, D, 1}, b::AbstractArray{T, D}, wb::W) where {D, T <: Real, W <: Union{Real, AbstractArray{<:Real, D}}}

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

@generated function _fit!(A::I, b::AbstractArray{T, D}, wb::W) where {P, D, T <: Real, I <: IndependentStatistic{T, D, P}, W <: Union{Real, AbstractArray{<:Real, D}}}
    code = Expr(:block)
    if W <: AbstractArray
        push!(
            code.args, quote
                size(b) == size(wb) || throw(ArgumentError("IndependentStatistic : size(b) != size(wb)"))
            end
        )
        WBI = :(wb[i])
    else
        push!(
            code.args, quote
                wb < 0 && throw(ArgumentError("weights can't be negative"))
                wb == 0 && return A
            end
        )
        WBI = :(wb)
    end
    push!(
        code.args, quote
            P > 1 || throw(ArgumentError("P must be greater than 1"))
            nonnegative(wb) || throw(ArgumentError("weights can't be negative"))
            wb == 0 && return A
            wa = weights(A)
            m1 = get_rawmoments(A, 1)
        end
    )
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
    push!(
        code.args, quote
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
        end
    )

    return code
end

function MomentsExpression(P::Int)
    code = Expr(:block)
    if P ≤ 2
        return code
    end
    for p in P:-1:3
        push!(code.args, :(get_rawmoments(A, $p)[i] += wai * BoN^$p + wbi * AoN^$p))
        for k in 1:(p - 2)
            push!(code.args, :(get_rawmoments(A, $p)[i] += binomial($p, $k) * get_rawmoments(A, $p - $k)[i] * BoN^$k))
        end
    end
    return code
end
