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

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse(string(
    "public STAT_HDU_KWD, STAT_GROUP_ID_KWD, STAT_MOMENT_INDEX_KWD, STAT_NB_MOMENTS_KWD, ",
    "STAT_WEIGHTS_KWD, isa_stat_hdu, find_stat_group_ids, find_stat_hdus")))

import OnlineStatsBase
import Statistics
import StatsBase

import OnlineStatsBase: _fit!, value

import StatsBase: fit!, kurtosis, mean, nobs, skewness, std, var, weights


include("UnivariateStatistics.jl")
include("IndependantStatistics.jl")
include("show.jl")

end
