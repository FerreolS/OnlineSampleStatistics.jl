using OnlineSampleStatistic
using Documenter

DocMeta.setdocmeta!(OnlineSampleStatistic, :DocTestSetup, :(using OnlineSampleStatistic); recursive=true)

makedocs(;
    modules=[OnlineSampleStatistic],
    authors="ferreol.soulez@univ-lyon1.fr",
    sitename="OnlineSampleStatistic.jl",
    format=Documenter.HTML(;
        canonical="https://FerreolS.github.io/OnlineSampleStatistic.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FerreolS/OnlineSampleStatistic.jl",
    devbranch="master",
)
