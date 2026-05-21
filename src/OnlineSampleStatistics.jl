"""
    OnlineSampleStatistics

Online, numerically stable sample statistics for scalar values and independent arrays.

The package provides mutable accumulators for first to higher-order moments with
single-pass updates, weighted updates, and merge operations.

This module returns exported constructors and statistics operations for online
accumulation.

# Example
```julia
using OnlineSampleStatistics
A = UnivariateStatistic(2)
fit!(A, [1.0, 2.0, 3.0])
mean(A)
```
"""
module OnlineSampleStatistics

export IndependentStatistic,
    UnivariateStatistic,
    fit!,
    get_moments,
    kurtosis,
    mean,
    merge!,
    nobs,
    order,
    skewness,
    std,
    var

VERSION >= v"1.11.0-DEV.469" && eval(
    Meta.parse(
        string(
            "public STAT_HDU_KWD, STAT_GROUP_ID_KWD, STAT_MOMENT_INDEX_KWD, STAT_NB_MOMENTS_KWD, ",
            "STAT_WEIGHTS_KWD, isa_stat_hdu, find_stat_group_ids, find_stat_hdus"
        )
    )
)

import OnlineStatsBase
import Statistics
import StatsBase

import OnlineStatsBase: _fit!

import StatsBase: fit!, kurtosis, mean, nobs, skewness, std, var, weights


include("UnivariateStatistics.jl")
include("IndependentStatistics.jl")
include("show.jl")

@doc """
    STAT_HDU_KWD

Keyword used to mark HDUs storing online statistics in FITS extensions.

This constant returns the FITS header keyword name as a `String`.
""" STAT_HDU_KWD

@doc """
    STAT_GROUP_ID_KWD

Keyword storing a group identifier for a set of statistics HDUs in FITS extensions.

This constant returns the FITS header keyword name as a `String`.
""" STAT_GROUP_ID_KWD

@doc """
    STAT_MOMENT_INDEX_KWD

Keyword storing the moment index in a statistics HDU in FITS extensions.

This constant returns the FITS header keyword name as a `String`.
""" STAT_MOMENT_INDEX_KWD

@doc """
    STAT_NB_MOMENTS_KWD

Keyword storing the number of tracked moments in a statistics HDU in FITS extensions.

This constant returns the FITS header keyword name as a `String`.
""" STAT_NB_MOMENTS_KWD

@doc """
    STAT_WEIGHTS_KWD

Keyword indicating the weights payload associated with a statistics HDU in FITS extensions.

This constant returns the FITS header keyword name as a `String`.
""" STAT_WEIGHTS_KWD

@doc """
    isa_stat_hdu(hdu)

Return whether a FITS HDU corresponds to serialized online statistics.

This function returns a `Bool`.
""" isa_stat_hdu

@doc """
    find_stat_group_ids(fitsfile)

Find available statistics group identifiers in a FITS file.

This function returns a collection of group identifiers.
""" find_stat_group_ids

@doc """
    find_stat_hdus(fitsfile, stat_group_id)

Find HDUs associated with a statistics group identifier in a FITS file.

This function returns the HDU tuple required to deserialize statistics.
""" find_stat_hdus

end
