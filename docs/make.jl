using DenoiseDNA
using Documenter

DocMeta.setdocmeta!(DenoiseDNA, :DocTestSetup, :(using DenoiseDNA); recursive=true)

makedocs(;
    modules=[DenoiseDNA],
    authors="Arthur Newbury",
    sitename="DenoiseDNA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/DenoiseDNA.jl",
        edit_link="master",
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/EvoArt/DenoiseDNA.jl",
    devbranch="master",
)
