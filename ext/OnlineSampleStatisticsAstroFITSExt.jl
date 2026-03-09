module OnlineSampleStatisticsAstroFITSExt

if isdefined(Base, :get_extension)
    using OnlineSampleStatistics, AstroFITS
    import OnlineSampleStatistics: isa_stat_hdu, find_stat_group_ids,
                                   STAT_HDU_KWD, STAT_GROUP_ID_KWD, STAT_MOMENT_INDEX_KWD,
                                   STAT_NB_MOMENTS_KWD, STAT_WEIGHTS_KWD
    import Base: read, write
else
    using ..OnlineSampleStatistics, ..AstroFITS
    import ..OnlineSampleStatistics: isa_stat_hdu, find_stat_group_ids,
                                     STAT_HDU_KWD, STAT_GROUP_ID_KWD, STAT_MOMENT_INDEX_KWD,
                                     STAT_NB_MOMENTS_KWD, STAT_WEIGHTS_KWD
    import ..Base: read, write
end

function isa_stat_hdu(hdu)
    (hdu isa AstroFITS.FitsImageHDU
     && haskey(hdu, STAT_HDU_KWD)
     && hdu[STAT_HDU_KWD].type == FITS_LOGICAL
     && hdu[STAT_HDU_KWD].logical)
end

function OnlineSampleStatistics.find_stat_group_ids(fitsfile::FitsFile)
    group_ids = Vector{String}()
    for hdu in fitsfile
        isa_stat_hdu(hdu) && push!(group_ids, hdu[STAT_GROUP_ID_KWD].string)
    end
    unique!(group_ids)
end

function Base.write(
    file::FitsFile,
    hdr::AstroFITS.OptionalHeader,
    stat::IndependentStatistic{T,N,K,W},
    stat_group_id ::String = string(rand('A':'Z', 16)...)
) where {T,N,K,W}
    dims = size(nobs(stat))

    moments_hdus = Vector{FitsImageHDU{T,N}}(undef, K)
    for k in 1:order(stat)
        moments_hdus[k] = FitsImageHDU{T,N}(file, dims)
        (k == 1) && merge!(moments_hdus[k], filter(!is_structural, hdr))
        merge!(moments_hdus[k], FitsHeader(
            STAT_HDU_KWD          => (true,          "is a OnlineSampleStatistics.jl data"),
            STAT_GROUP_ID_KWD     => (stat_group_id, "ID to group statistical moments HDUs"),
            STAT_NB_MOMENTS_KWD   => (K,             "number of statistical moments"),
            STAT_MOMENT_INDEX_KWD => (k,             "th statistical moment")))
        # adding EXTNAME unless the user already specified one
        haskey(moments_hdus[k], "EXTNAME") || push!(moments_hdus[k], "EXTNAME" => "MOMENT-$k")
        write(moments_hdus[k], OnlineSampleStatistics.get_rawmoments(stat, k))
    end
    
    weights_hdu = FitsImageHDU{W,N}(file, dims)
    merge!(weights_hdu, FitsHeader(
        STAT_HDU_KWD      => (true,          "is a OnlineSampleStatistics.jl data"),
        STAT_GROUP_ID_KWD => (stat_group_id, "ID to group statistical moments HDUs"),
        STAT_WEIGHTS_KWD  => (true,          "is statistical weights")))
    # adding EXTNAME unless the user already specified one
    haskey(weights_hdu, "EXTNAME") || push!(weights_hdu, "EXTNAME" => "WEIGHTS")
    write(weights_hdu, nobs(stat))
    
    return file
end

function find_stat_hdus(fitsfile::FitsFile, stat_group_id::String)
    local moments_hdus, weights_hdu, T, N, K, W
    for hdu in fitsfile
        isa_stat_hdu(hdu) || continue
        hdu[STAT_GROUP_ID_KWD].string == stat_group_id || continue

        if haskey(hdu, STAT_WEIGHTS_KWD) && hdu[STAT_WEIGHTS_KWD].logical
            weights_hdu = hdu
            W = weights_hdu.data_eltype
        else
             if !(@isdefined moments_hdus)
                T = hdu.data_eltype
                N = hdu.data_ndims
                K = hdu[STAT_NB_MOMENTS_KWD].integer
                moments_hdus = Vector{FitsImageHDU}(undef, K)
            end
            T = promote_type(T, hdu.data_eltype)
            k = hdu[STAT_MOMENT_INDEX_KWD].integer
            moments_hdus[k] = hdu
        end
    end
    (@isdefined weights_hdu)  || throw(ArgumentError("could not find weights HDU"))
    (@isdefined moments_hdus) || throw(ArgumentError("could not find any moment HDU"))
    for k in 1:K
        isassigned(moments_hdus, k) || throw(ArgumentError("could not find moment number $k HDU"))
    end
    (moments_hdus, weights_hdu, T, N, K, W)
end

function Base.read(
    ::Type{IndependentStatistic},
    fitsfile::FitsFile,
    stat_group_id::String
    ; readkwds...
)
    (moments_hdus, weights_hdu, T, N, K, W) = find_stat_hdus(fitsfile, stat_group_id)

    moments = NTuple{K,Array{T,N}}( read(Array{T,N}, moments_hdus[k]) for k in 1:K )
    weights = read(Array{W,N}, weights_hdu; readkwds...)

    OnlineSampleStatistics.build_from_rawmoments(weights, moments)
end

function Base.read(
    ::Type{IndependentStatistic},
    fitsfile::FitsFile
    ; ext::Union{Int,String}=1,
      readkwds...)
    isa_stat_hdu(fitsfile[ext]) || throw(ArgumentError("HDU \"$ext\" is not a stat HDU"))
    stat_group_id = fitsfile[ext][STAT_GROUP_ID_KWD].string
    read(IndependentStatistic, fitsfile, stat_group_id; readkwds...)
end

end