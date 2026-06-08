import CompariMotif
import Documenter
import DocumenterCitations

Documenter.DocMeta.setdocmeta!(
    CompariMotif,
    :DocTestSetup,
    :(import CompariMotif);
    recursive = true
)
bib = DocumenterCitations.CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    style = :authoryear
)

Documenter.makedocs(;
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
        "Regex Syntax" => "regex_syntax.md",
        "FAQ / How-To" => "faq.md",
        "Internal API & Pipeline" => "internal_api.md",
        "References" => "references.md"
    ]
)

Documenter.deploydocs(;
    repo = "github.com/diegozea/CompariMotif.jl.git",
    devbranch = "main"
)
