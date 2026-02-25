using CompariMotif
using Documenter

DocMeta.setdocmeta!(CompariMotif, :DocTestSetup, :(using CompariMotif); recursive=true)

makedocs(;
    modules=[CompariMotif],
    authors="Diego Javier Zea <diegozea@gmail.com> and contributors",
    sitename="CompariMotif.jl",
    format=Documenter.HTML(;
        canonical="https://diegozea.github.io/CompariMotif.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/diegozea/CompariMotif.jl",
    devbranch="main",
)
