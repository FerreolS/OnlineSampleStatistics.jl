module OnlineSampleStatistics
using ZippedArrays

export UnivariateStatistic,
    mean,
    var,
    nobs,
    skewness,
    kurtosis,
    weights


include("UnivariateStat.jl")

end
