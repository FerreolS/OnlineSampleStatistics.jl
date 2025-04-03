module OnlineSampleStatistics

export UnivariateStatistic,
    mean,
    var,
    nobs,
    skewness,
    kurtosis,
    weights,
    IndependentStatistic


include("UnivariateStat.jl")
include("IndependantStatistics.jl")

end
