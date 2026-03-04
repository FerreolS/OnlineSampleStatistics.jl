module OnlineSampleStatisticsAstroFITSExt

if isdefined(Base, :get_extension)
    using OnlineSampleStatistics, AstroFITS, FITSHeaders
    import OnlineSampleStatistics: STAT_HDU_KWD, isa_stat_hdu, weights
    import AstroFITS: FitsImageHDU, OptionalHeader, type_to_bitpix
    import Base: read, write
else
    using ..OnlineSampleStatistics, ..AstroFITS, ..FITSHeaders
    import ..OnlineSampleStatistics: STAT_HDU_KWD, isa_stat_hdu, weights
    import ..AstroFITS: FitsImageHDU, OptionalHeader, type_to_bitpix
    import ..Base: read, write
end

function __init__()
    println(stderr, methods(AstroFITS.type_to_bitpix))
end

function OnlineSampleStatistics.isa_stat_hdu(hdu)
    (hdu isa AstroFITS.FitsImageHDU
     && haskey(hdu, STAT_HDU_KWD)
     && hdu[STAT_HDU_KWD].type == FITSHeaders.FITS_LOGICAL
     && hdu[STAT_HDU_KWD].logical)
end



#"add a stat image HDU to a FITS file, then add other HDUs"
#function Base.write(file::FitsFile, header::OptionalHeader, data::IndependentStatistic, args...)
#    write(file, header, data)
#    write(file, args...)
#end
#
#"add a stat image HDU to a FITS file"
#function Base.write(file::FitsFile, hdr::OptionalHeader, stat::IndependentStatistic)
#    write(merge!(FitsImageHDU(file, stat), hdr), stat)
#    return file # returns the file not the HDU
#end
#
#
#
#"""
#    FitsImageHDU{T,N}(::FitsFile, ::IndependentStatistic{<:Real}) -> FitsImageHDU
#
#declare a stat image HDU in a FITS file, from an IndependentStatistic
#
#last dimension is for the moments, with weights at the end.
#"""
#function AstroFITS.FitsImageHDU{T,Np1}(file::FitsFile,
#                                       stat::IndependentStatistic{<:Real,N}
#) where {T<:Real,Np1,N}
#    Np1 == N+1 || throw(DimensionMismatch("HDU has $Np1 dimensions instead of $(N+1)")
#    dims = (size(weights(stat))..., order(stat)+1)  # +1 for the weights
#    return FitsImageHDU{T,Np1}(file, dims)
#end
#
#"alias with only type {T} specified"
#function AstroFITS.FitsImageHDU{T}(
#    file::FitsFile, stat::IndependentStatistic{<:Real,N}
#) where {T<:Real,N}
#    FitsImageHDU{T,N+1}(file, stat) # +1 dimension for the moments
#end
#
#"alias with no type specified"
#function AstroFITS.FitsImageHDU(
#    file::FitsFile, stat::IndependentStatistic{T,N,K,W}
#) where {T<:Real,N,K,W<:Real}
#    FitsImageHDU{promote_type(T,W)}(file, stat)
#end

function OnlineSampleStatistics.declare_HDUs(
    file::FitsFile, stat::IndependentStatistic{T,N,K,W}) where {T,N,K,W}
    println(stderr, methods(AstroFITS.type_to_bitpix))
    moments_hdus = Vector{FitsImageHDU{T,N}}(undef, K)
    dims = size(weights(stat))
    for k in 1:order(stat)
        moments_hdus[k] = FitsImageHDU{T,N}(file, dims)
    end
    weights_hdu = FitsImageHDU{K,N}(file, dims)
    (moments_hdus, weights_hdu)
end

const STAT_GROUP_ID_KWD = "STAT-GROUP-ID"
const STAT_MOMENT_INDEX_KWD = "STAT-MOMENT-INDEX"
const STAT_NB_MOMENTS_KWD = "STAT-NB-MOMENTS"
const STAT_WEIGHTS_KWD = "STAT-WEIGHTS"

"""
    write(moments_hdus::Vector{FitsImageHDU}, weights_hdu::FitsImageHDU, stat::IndependentStatistic) -> (moments_hdus, weights_hdu)

write IndependentStatistic' raw moments and weights to images HDUs 
"""
function Base.write(
    moments_hdus::Vector{FitsImageHDU}, weights_hdu::FitsImageHDU, stat::IndependentStatistic)

    K = order(stat)
    length(moments_hdu) == K || throw(DimensionMismatch(
        "Number of moments HDUs is $(length(moments_hdus)) instead of $K"))

    # ID to be able to retrieve the different HDUs inside the FITS file, even if they get reordered
    random_id = round(Int, randn() * 10^12)

    for k in 1:K
        push!(moments_hdu[k],
            STAT_HDU_KWD          => (true,      "is a OnlineSampleStatistics.jl data"),
            STAT_GROUP_ID_KWD     => (random_id, "stat ID to group moments HDUs together"),
            STAT_NB_MOMENTS_KWD   => (K,         "number of IndependentStatistic moments"),
            STAT_MOMENT_INDEX_KWD => (k,         "th IndependentStatistic moment"))
        write(moments_hdu[k], get_rawmoments(stat, k))
    end

    push!(weights_hdu,
        STAT_HDU_KWD        => (true,      "OnlineSampleStatistics.jl data"),
        STAT_GROUP_ID_KWD   => (random_id, "stat ID to group moments HDUs together"),
        STAT_WEIGHTS_KWD    => (true,      "IndependentStatistic weights"))
    write(weights_hdu, weights(stat))

    return (moments_hdus, weights_hdu)
end

function Base.read(
    ::Type{IndependentStatistic}, fitsfile::FitsFile, stat_group_id::Int
    ; readkwds...)

    local T, N, K, moments, weights # will be filled when reading HDUs

    for hdu in fitsfile
        # read only stats HDUs with the given group id
        isa_stat_hdu(hdu) || continue
        hdu[STAT_GROUP_ID_KWD].integer == stat_group_id || continue

        if haskey(hdu, STAT_WEIGHTS_KWD) && hdu[STAT_WEIGHTS_KWD].logical
            if !(@isdefined moments)
                T = hdu.data_eltype
                N = hdu.data_ndims
                K = hdu[STAT_NB_MOMENTS_KWD].integer
                moments = Vector{Array{T,N}}(undef, K)
            end
            k = hdu[STAT_MOMENT_KWD].integer
            moments[k] = read(Array{T,N}, hdu; readkwds...)
        else
            weights = read(hdu; readkwds...)
        end

        # stop reading HDUs if work is done
        if (@isdefined weights) && (@isdefined moments) && all(i -> isassigned(moments, i), 1:K)
            break
        end
    end
    
    # create an IndependentStatistic from raw moments
    OnlineSampleStatistics.build_from_rawmoments(weights, moments)
end

function Base.read(
    ::Type{IndependentStatistic}, fitsfile::FitsFile, hdu_ext::Union{Int,String}=1
    ; readkwds...)
    stat_group_id = fitsfile[hdu_ext][STAT_GROUP_ID_KWD].integer
    Base.read(IndependentStatistic, fitsfile, stat_group_id)
end

    
#
#"""
#    write(::FitsImageHDU, stat::IndependentStatistic{<:Real}) -> FitsImageHDU
#
#write IndependentStatistic' raw moments and weights to a stat image HDU
#"""
#function Base.write(
#    hdu::FitsImageHDU{T,Np1}, stat::IndependentStatistic{T2,N,K,W}
#) where {T<:Real,Np1,T2<:Real,N,K<:Real}
#
#    Np1 == N+1 || throw(DimensionMismatch("HDU has $Np1 dimensions instead of $(N+1)"))
#    dims = (size(weights(stat))..., order(stat)+1)  # +1 for the weights
#    hdu.data_size == dims || throw(
#        DimensionMismatch("HDU size is $(hdu.data_size) instead of $dims"))
#
#    bitpix_weights = AstroFITS.type_to_bitpix(W)
#
#    push!(hdu, STAT_HDU_KWD => (true, "image is IndependentStatistic data"),
#               BITPIX_WEIGHTS_KWD => (bitpix_weights, "initial BITPIX for weights")
#               BITPIX_WEIGHTS_KWD => (bitpix_weights, "initial BITPIX for weights"))
#
#    i = 1
#    for k in 1:K
#        momentk = get_rawmoments(stat, k)
#        write(hdu, momentk; first=i)
#        i += length(momentk)
#    end
#    write(hdu, weights(args); first=i)
#
#    return hdu
#end
#
#
#"""
#    read(::Type{IndependentStatistic{T,N,K,W}, hdu::FitsImageHDU; readkwds...) -> IndependentStatistic
#
#read IndependentStatistic' raw moments and weights from a stat image HDU
#"""
#function Base.read(::Type{IndependentStatistic{T,N,K,W},
#                   hdu::FitsImageHDU{T2,Np1}
#                   ; readkwds...
#) where {T,N,K,W,T2,Np1}
#
#    Np1 == N+1 || throw(DimensionMismatch("HDU has $Np1 dimensions instead of $(N+1)"))
#
#    K == hdu.data_size[end]-1 || throw(
#        DimensionMismatch("HDU has $(hdu.data_size[end]-1) statistical moments instead of $K")
#
#    n = hdu[STAT_NB_SAMPLES_KWD].integer
#    data = read(Array{T,N1}, hdu; readkwds...)
#    s = NTuple{L,A}( A(sl) for sl in eachslice(data; dims=N1))
#
#    return build_from_rawmoments  IndependentStatistic{L,T,N,A}(s,n)
#end
#
#
#Base.read(::Type{IndependentStatistic}, hdu::FitsImageHDU; kwds...) =
#    read(IndependentStatistic{hdu.data_size[end]}, hdu; kwds...)
#
#Base.read(::Type{IndependentStatistic{L}}, hdu::FitsImageHDU{T}; kwds...) where {L,T} =
#    read(IndependentStatistic{L,T}, hdu; kwds...)
#
#Base.read(::Type{IndependentStatistic{L,T}}, hdu::FitsImageHDU{<:Any,N1}; kwds...) where {L,T,N1} =
#    read(IndependentStatistic{L,T,N1-1}, hdu; kwds...)
#
#Base.read(::Type{IndependentStatistic{L,T,N}}, hdu::FitsImageHDU; kwds...) where {L,T,N} =
#    read(IndependentStatistic{L,T,N,Array{T,N}}, hdu; kwds...)


end