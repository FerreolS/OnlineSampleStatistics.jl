using OnlineSampleStatistics
using AstroFITS
using Documenter

const OnlineSampleStatisticsAstroFITSExt = Base.get_extension(
    OnlineSampleStatistics,
    :OnlineSampleStatisticsAstroFITSExt
)
@assert !isnothing(OnlineSampleStatisticsAstroFITSExt)
@eval OnlineSampleStatistics const OnlineSampleStatisticsAstroFITSExt = $OnlineSampleStatisticsAstroFITSExt

DocMeta.setdocmeta!(OnlineSampleStatistics, :DocTestSetup, :(using OnlineSampleStatistics); recursive = true)

makedocs(;
    modules = [OnlineSampleStatistics, OnlineSampleStatisticsAstroFITSExt],
    authors = "ferreol.soulez@univ-lyon1.fr",
    sitename = "OnlineSampleStatistics.jl",
    format = Documenter.HTML(;
        canonical = "https://FerreolS.github.io/OnlineSampleStatistics.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Core Types" => "types.md",
    ],
)

deploydocs(;
    repo = "github.com/FerreolS/OnlineSampleStatistics.jl",
    devbranch = "master",
)
