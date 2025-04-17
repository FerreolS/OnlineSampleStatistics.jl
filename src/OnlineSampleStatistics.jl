module OnlineSampleStatistics

export UnivariateStatistic,
    IndependentStatistic,
    mean,
    var,
    std,
    nobs,
    skewness,
    kurtosis,
    weights,
    order,
    get_moments

import OnlineStatsBase: value, fit!



include("UnivariateStatistics.jl")
include("IndependantStatistics.jl")

end
