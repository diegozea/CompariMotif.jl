using TestItems

@testitem "external API docs examples" begin
    using Test
    using CompariMotif

    # This test mirrors the examples in `docs/src/external_api.md`.
    # Keep this file and that page in sync when changing the example inputs,
    # the documented outputs, or the public API behavior.

    @testset "quick start" begin
        options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
        result = compare("RKLI", "R[KR]L[IV]", options)

        @test (
            result.matched,
            result.query_relationship,
            result.search_relationship,
            result.matched_positions
        ) == (true, "Variant Match", "Degenerate Match", 4)
    end

    @testset "matrix comparisons" begin
        motifs = ["RKLI", "R[KR]L[IV]", "RxLE"]
        options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
        results = compare(motifs, options)
        table = to_column_table(results)

        @test size(results) == (3, 3)
        @test length(table.query) == 9
    end

    @testset "custom residue frequencies" begin
        weighted = ComparisonOptions(;
            alphabet = DNAAlphabet(),
            residue_frequencies = Dict('A' => 0.3, 'C' => 0.2, 'G' => 0.2, 'T' => 0.3),
            min_shared_positions = 1,
            normalized_ic_cutoff = 0.0
        )
        result = compare("ATG", "[AGT]TG", weighted)

        @test round(result.match_ic, digits = 3) == 2.19
    end

    @testset "canonicalization" begin
        @test normalize_motif("r[kR].{0,1}l") == "R[RK]x{0,1}L"
    end
end
