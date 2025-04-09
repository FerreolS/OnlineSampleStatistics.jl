module OnlineSampleStatistics

export UnivariateStatistic,
    mean,
    var,
    nobs,
    skewness,
    kurtosis,
    weights,
    order,
    IndependentStatistic


include("UnivariateStat.jl")
include("IndependantStatistics.jl")

end
