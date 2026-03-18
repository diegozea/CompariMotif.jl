@testitem "paper Figure 2 examples" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    # Edwards et al. (2008), Bioinformatics 24(10):1307-1309, Fig. 2 (page 1308)
    query = "[KR]xLx[FYLIMVP]"
    cases = [
        ("[KR]xLx[FYLIMVP]", "Exact Match", "Exact Match"),
        ("[KR]xL", "Exact Parent", "Exact Subsequence"),
        ("x[KR]xL", "Exact Overlap", "Exact Overlap"),
        ("RxLPL", "Degenerate Match", "Variant Match"),
        ("RxLE", "Degenerate Parent", "Variant Subsequence"),
        ("[KR]x[IL]", "Variant Parent", "Degenerate Subsequence"),
        ("LxPP", "Degenerate Overlap", "Variant Overlap"),
        ("RSxPP", "Complex Match", "Complex Match"),
        ("L[IL]xL", "Complex Parent", "Complex Subsequence"),
        ("RxRS[ILMV]", "Complex Overlap", "Complex Overlap")
    ]

    for (search, query_relationship, search_relationship) in cases
        @testset "$query vs $search" begin
            result = compare(query, search, options)
            @test result.matched
            @test result.query_relationship == query_relationship
            @test result.search_relationship == search_relationship
        end
    end
end
