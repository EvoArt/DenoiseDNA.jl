using DenoiseDNA
using Documenter

DocMeta.setdocmeta!(DenoiseDNA, :DocTestSetup, :(using DenoiseDNA); recursive=true)

makedocs(;
    modules=[DenoiseDNA],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/DenoiseDNA.jl/blob/{commit}{path}#{line}",
    sitename="DenoiseDNA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/DenoiseDNA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/DenoiseDNA.jl",
)
