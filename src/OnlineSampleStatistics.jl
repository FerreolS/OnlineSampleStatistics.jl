module OnlineSampleStatistics

export UnivariateStatistic,
    IndependentStatistic,
    mean,
    var,
    std,
    nobs,
    skewness,
    kurtosis,
    order,
    get_moments,
    fit!,
    merge!

VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse(string(
    "public STAT_HDU_KWD, STAT_GROUP_ID_KWD, STAT_MOMENT_INDEX_KWD, STAT_NB_MOMENTS_KWD, ",
    "STAT_WEIGHTS_KWD, isa_stat_hdu, find_stat_group_ids, find_stat_hdus")))

import OnlineStatsBase
import Statistics
import StatsBase

import OnlineStatsBase: value, _fit!

import StatsBase: fit!, nobs, mean, var, std, skewness, kurtosis, weights


include("UnivariateStatistics.jl")
include("IndependantStatistics.jl")
include("show.jl")

end
