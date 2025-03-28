using OnlineSampleStatistics
using Documenter

DocMeta.setdocmeta!(OnlineSampleStatistics, :DocTestSetup, :(using OnlineSampleStatistics); recursive=true)

makedocs(;
    modules=[OnlineSampleStatistics],
    authors="ferreol.soulez@univ-lyon1.fr",
    sitename="OnlineSampleStatistics.jl",
    format=Documenter.HTML(;
        canonical="https://FerreolS.github.io/OnlineSampleStatistics.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FerreolS/OnlineSampleStatistics.jl",
    devbranch="master",
)
