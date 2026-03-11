using TestItems

@testitem "paper Table 1 relationship words" begin
    using Test
    using CompariMotif

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    # Edwards et al. (2008), Table 1 (page 1308)
    cases = [
        ("RKLI", "RKLI", "Exact Match"),
        ("ARKLI", "RKLI", "Exact Parent"),
        ("RKLI", "ARKLI", "Exact Subsequence"),
        ("ARKLI", "RKLIT", "Exact Overlap"),
        ("RKLI", "R[KR]L[IV]", "Variant Match"),
        ("ARKLI", "R[KR]L[IV]", "Variant Parent"),
        ("RKLI", "AR[KR]L[IV]", "Variant Subsequence"),
        ("ARKLI", "R[KR]L[IV]T", "Variant Overlap"),
        ("R[KR]L[IV]", "RKLI", "Degenerate Match"),
        ("AR[KR]L[IV]", "RKLI", "Degenerate Parent"),
        ("R[KR]L[IV]", "ARKLI", "Degenerate Subsequence"),
        ("AR[KR]L[IV]", "RKLIT", "Degenerate Overlap"),
        ("[KR]L[IMV]", "[RKQ]L[IV]", "Complex Match"),
        ("A[KR]L[IMV]", "[RKQ]L[IV]", "Complex Parent"),
        ("[KR]L[IMV]", "A[RKQ]L[IV]", "Complex Subsequence"),
        ("A[KR]L[IMV]", "[RKQ]L[IV]T", "Complex Overlap")
    ]

    for (query, search, relationship) in cases
        result = compare(query, search, options)
        @test result.matched
        @test result.query_relationship == relationship
    end
end
