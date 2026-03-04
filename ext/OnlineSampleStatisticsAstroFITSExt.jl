module OnlineSampleStatisticsAstroFITSExt

if isdefined(Base, :get_extension)
    using OnlineSampleStatistics, AstroFITS, FITSHeaders
    import OnlineSampleStatistics: STAT_HDU_KWD, weights, find_stat_groupd_ids
    import Base: read, write
else
    using ..OnlineSampleStatistics, ..AstroFITS, ..FITSHeaders
    import ..OnlineSampleStatistics: STAT_HDU_KWD, weights, find_stat_groupd_ids
    import ..Base: read, write
end

function isa_stat_hdu(hdu)
    (hdu isa AstroFITS.FitsImageHDU
     && haskey(hdu, STAT_HDU_KWD)
     && hdu[STAT_HDU_KWD].type == FITSHeaders.FITS_LOGICAL
     && hdu[STAT_HDU_KWD].logical)
end

function declare_HDUs(
    file::FitsFile, stat::IndependentStatistic{T,N,K,W}) where {T,N,K,W}
    moments_hdus = Vector{FitsImageHDU{T,N}}(undef, K)
    dims = size(weights(stat))
    for k in 1:order(stat)
        moments_hdus[k] = FitsImageHDU{T,N}(file, dims)
    end
    weights_hdu = FitsImageHDU{W,N}(file, dims)
    (moments_hdus, weights_hdu)
end

const STAT_GROUP_ID_KWD = "STAT-GROUP-ID"
const STAT_MOMENT_INDEX_KWD = "STAT-MOMENT-INDEX"
const STAT_NB_MOMENTS_KWD = "STAT-NB-MOMENTS"
const STAT_WEIGHTS_KWD = "STAT-WEIGHTS"

function OnlineSampleStatistics.find_stat_groupd_ids(fitsfile::FitsFile)
    group_ids = Set{String}()
    for hdu in fitsfile
        isa_stat_hdu(hdu) && push!(group_ids, hdu[STAT_GROUP_ID_KWD].string)
    end
    group_ids
end

function Base.write(
    file::FitsFile,
    hdr::AstroFITS.OptionalHeader,
    stat::IndependentStatistic,
    stat_group_id ::String = string(rand('A':'Z', 16)...)
)
    (moments_hdus, weights_hdu) = declare_HDUs(file, stat)
    merge!(moments_hdus[1], filter(!is_structural, hdr))
    write(moments_hdus, weights_hdu, stat, stat_group_id)
    return file # returns the file not the HDU
end

"""
    write(moments_hdus::Vector{FitsImageHDU{T,N}}, weights_hdu::FitsImageHDU{W,N}, stat::IndependentStatistic; stat_group_id::Int=rand()) -> (moments_hdus, weights_hdu)

write IndependentStatistic' raw moments and weights to images HDUs

`stat_group_id`: ID to be able to retrieve the different HDUs inside the FITS file, even if they get reordered
"""
function Base.write(
    moments_hdus::Vector{FitsImageHDU{T,N}},
    weights_hdu::FitsImageHDU{W,N},
    stat::IndependentStatistic,
    stat_group_id::String,
) where {T,N,W}

    K = order(stat)
    length(moments_hdus) == K || throw(DimensionMismatch(
        "Number of moments HDUs is $(length(moments_hdus)) instead of $K"))

    for k in 1:K
        push!(moments_hdus[k],
            STAT_HDU_KWD          => (true,          "is a OnlineSampleStatistics.jl data"),
            STAT_GROUP_ID_KWD     => (stat_group_id, "stat ID to group moments HDUs together"),
            STAT_NB_MOMENTS_KWD   => (K,             "number of IndependentStatistic moments"),
            STAT_MOMENT_INDEX_KWD => (k,             "th IndependentStatistic moment"))
        write(moments_hdus[k], OnlineSampleStatistics.get_rawmoments(stat, k))
    end

    push!(weights_hdu,
        STAT_HDU_KWD      => (true,          "OnlineSampleStatistics.jl data"),
        STAT_GROUP_ID_KWD => (stat_group_id, "stat ID to group moments HDUs together"),
        STAT_WEIGHTS_KWD  => (true,          "IndependentStatistic weights"))
    write(weights_hdu, weights(stat))

    return (moments_hdus, weights_hdu)
end

function Base.read(
    ::Type{IndependentStatistic}, fitsfile::FitsFile, stat_group_id::String
    ; readkwds...)

    local T, N, K, moments, weights # will be filled when reading HDUs

    for hdu in fitsfile
        # read only stats HDUs with the given group id
        isa_stat_hdu(hdu) || continue
        hdu[STAT_GROUP_ID_KWD].string == stat_group_id || continue

        if haskey(hdu, STAT_WEIGHTS_KWD) && hdu[STAT_WEIGHTS_KWD].logical
            weights = read(hdu; readkwds...)
        else
            if !(@isdefined moments)
                T = hdu.data_eltype
                N = hdu.data_ndims
                K = hdu[STAT_NB_MOMENTS_KWD].integer
                moments = Vector{Array{T,N}}(undef, K)
            end
            k = hdu[STAT_MOMENT_INDEX_KWD].integer
            moments[k] = read(Array{T,N}, hdu; readkwds...)
        end

        # stop reading HDUs if work is done
        if (@isdefined weights) && (@isdefined moments) && all(i -> isassigned(moments, i), 1:K)
            break
        end
    end

    if !(@isdefined weights)
        error("could not find weights HDU")
    end
    for k in 1:K
        if !(@isdefined moments) || !isassigned(moments, k)
            error("could not find moment number $k HDU")
        end
    end

    moments = tuple(moments...) # MANDATORY

    # create an IndependentStatistic from raw moments
    OnlineSampleStatistics.build_from_rawmoments(weights, moments)
end

function Base.read(
    ::Type{IndependentStatistic}, fitsfile::FitsFile; ext::Union{Int,String}=1, readkwds...)
    stat_group_id = fitsfile[ext][STAT_GROUP_ID_KWD].string
    Base.read(IndependentStatistic, fitsfile, stat_group_id; readkwds...)
end

end