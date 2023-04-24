using Erosion
using Documenter

DocMeta.setdocmeta!(Erosion, :DocTestSetup, :(using Erosion); recursive=true)

makedocs(;
    modules=[Erosion],
    authors="Cédric BELMANT",
    repo="https://github.com/serenity4/Erosion.jl/blob/{commit}{path}#{line}",
    sitename="Erosion.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://serenity4.github.io/Erosion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/serenity4/Erosion.jl",
    devbranch="main",
)
