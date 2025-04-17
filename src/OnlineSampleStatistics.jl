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
    get_moments

import OnlineStatsBase: value, _fit!, _merge!

import StatsBase: fit!, nobs, mean, var, std, skewness, kurtosis, weights



include("UnivariateStatistics.jl")
include("IndependantStatistics.jl")

end
