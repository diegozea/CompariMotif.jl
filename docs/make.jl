using CompariMotif
using Documenter

DocMeta.setdocmeta!(CompariMotif, :DocTestSetup, :(using CompariMotif); recursive = true)

makedocs(;
    modules = [CompariMotif],
    checkdocs = :exports,
    authors = "Diego Javier Zea <diegozea@gmail.com> and contributors",
    sitename = "CompariMotif.jl",
    format = Documenter.HTML(;
        canonical = "https://diegozea.github.io/CompariMotif.jl",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "External API" => "external_api.md",
        "FAQ / How-To" => "faq.md",
        "Internal API & Pipeline" => "internal_api.md"
    ]
)

deploydocs(;
    repo = "github.com/diegozea/CompariMotif.jl",
    devbranch = "main"
)
