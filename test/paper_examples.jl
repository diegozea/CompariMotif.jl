using TestItems

@testitem "paper Figure 1 example" begin
    using Test
    using CompariMotif

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    # Edwards et al. (2008), Bioinformatics 24(10):1307-1309, Fig. 1 (page 1308)
    result = compare("[KR]xLx[FYLIMVP]", "RxLE", options)

    @test result.matched
    @test result.query_relationship == "Degenerate Parent"
    @test result.search_relationship == "Variant Subsequence"
    @test result.matched_positions == 2
    @test result.match_ic ≈ 1.769 atol = 1e-3
    @test result.normalized_ic ≈ 0.835 atol = 1e-3
    @test result.core_ic ≈ 0.590 atol = 1e-3
    @test result.score ≈ 1.669 atol = 1e-3
end

@testitem "paper Table 1 relationship codes" begin
    using Test
    using CompariMotif

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    # Edwards et al. (2008), Table 1 (page 1308)
    cases = [
        ("RKLI", "RKLI", "e-m"),
        ("ARKLI", "RKLI", "e-p"),
        ("RKLI", "ARKLI", "e-s"),
        ("ARKLI", "RKLIT", "e-o"),
        ("RKLI", "R[KR]L[IV]", "v-m"),
        ("ARKLI", "R[KR]L[IV]", "v-p"),
        ("RKLI", "AR[KR]L[IV]", "v-s"),
        ("ARKLI", "R[KR]L[IV]T", "v-o"),
        ("R[KR]L[IV]", "RKLI", "d-m"),
        ("AR[KR]L[IV]", "RKLI", "d-p"),
        ("R[KR]L[IV]", "ARKLI", "d-s"),
        ("AR[KR]L[IV]", "RKLIT", "d-o"),
        ("[KR]L[IMV]", "[RKQ]L[IV]", "c-m"),
        ("A[KR]L[IMV]", "[RKQ]L[IV]", "c-p"),
        ("[KR]L[IMV]", "A[RKQ]L[IV]", "c-s"),
        ("A[KR]L[IMV]", "[RKQ]L[IV]T", "c-o")
    ]

    for (query, search, code) in cases
        result = compare(query, search, options)
        @test result.matched
        @test result.query_relationship_code == code
    end
end
