@testitem "documentation doctests" begin
    import Documenter

    Core.eval(Main, :(import CompariMotif))
    Documenter.DocMeta.setdocmeta!(
        CompariMotif,
        :DocTestSetup,
        :(using CompariMotif);
        recursive = true
    )

    Documenter.doctest(CompariMotif; manual = true)
end
