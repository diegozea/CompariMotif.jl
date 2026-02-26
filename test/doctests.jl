using TestItems

@testitem "documentation doctests" begin
    using Documenter
    using CompariMotif

    Documenter.doctest(CompariMotif)
end
