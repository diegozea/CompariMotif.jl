using TestItems

@testitem "documentation doctests" begin
    using Documenter
    using CompariMotif

    Core.eval(Main, :(import CompariMotif))
    DocMeta.setdocmeta!(CompariMotif, :DocTestSetup, :(using CompariMotif); recursive = true)

    Documenter.doctest(CompariMotif; manual = true)
end
