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

import OnlineStatsBase
import Statistics
import StatsBase

import OnlineStatsBase: value, _fit!

import StatsBase: fit!, nobs, mean, var, std, skewness, kurtosis, weights


include("UnivariateStatistics.jl")
include("IndependantStatistics.jl")
include("show.jl")

end
