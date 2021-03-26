using PyThermo
using Documenter

DocMeta.setdocmeta!(PyThermo, :DocTestSetup, :(using PyThermo); recursive=true)

makedocs(;
    modules=[PyThermo],
    authors="Alex Ames <alexander.m.ames@gmail.com> and contributors",
    repo="https://github.com/stillyslalom/PyThermo.jl/blob/{commit}{path}#{line}",
    sitename="PyThermo.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stillyslalom.github.io/PyThermo.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stillyslalom/PyThermo.jl",
)
