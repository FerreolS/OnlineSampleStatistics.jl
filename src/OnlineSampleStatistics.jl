module OnlineSampleStatistics

export UnivariateStatistic,
    IndependentStatistic,
    mean,
    var,
    nobs,
    skewness,
    kurtosis,
    weights,
    order


include("UnivariateStatistics.jl")
include("IndependantStatistics.jl")

end
