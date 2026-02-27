using TestItems

@testitem "README minimal example" begin
    using Test
    using CompariMotif

    motifs = ["RKLI", "R[KR]L[IV]", "[KR]xLx[FYLIMVP]", "RxLE"]
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    results = compare(motifs, options)

    @test size(results) == (4, 4)
    @test results[3, 4].matched
    @test results[3, 4].query_relationship == "Degenerate Parent"
    @test results[3, 4].matched_positions == 2

    table = to_column_table(results)
    @test length(table.query_index) == 16
    @test :query_relationship in keys(table)
end
