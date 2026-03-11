using CompariMotif
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(CompariMotif, :DocTestSetup, :(using CompariMotif); recursive = true)
bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :authoryear)

makedocs(;
    modules = [CompariMotif],
    checkdocs = :exports,
    authors = "Diego Javier Zea <diegozea@gmail.com> and contributors",
    sitename = "CompariMotif.jl",
    format = Documenter.HTML(;
        canonical = "https://diegozea.github.io/CompariMotif.jl",
        edit_link = "main",
        assets = String["assets/citations.css"]
    ),
    plugins = [bib],
    pages = [
        "Home" => "index.md",
        "External API" => "external_api.md",
        "FAQ / How-To" => "faq.md",
        "Internal API & Pipeline" => "internal_api.md",
        "References" => "references.md"
    ]
)

deploydocs(;
    repo = "github.com/diegozea/CompariMotif.jl",
    devbranch = "main"
)
