using TestItems

@testitem "README minimal example" begin
    using Test
    using CompariMotif

    motifs = ["RKLI", "R[KR]L[IV]", "[KR]xLx[FYLIMVP]", "RxLE"]
    results = compare(motifs; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    @test size(results) == (4, 4)
    @test results[3, 4].matched
    @test results[3, 4].query_relationship == "Degenerate Parent"
    @test results[3, 4].matched_positions == 2

    outfile = joinpath(mktempdir(), "comparimotif_results.tsv")
    write_results_tsv(outfile, motifs, motifs, results)
    @test isfile(outfile)
end
