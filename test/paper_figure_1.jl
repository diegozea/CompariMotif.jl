@testitem "paper Figure 1 example" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    # Edwards et al. (2008), Bioinformatics 24(10):1307-1309, Fig. 1.
    # The numeric checks below are the figure's printed values.
    result = compare("[KR]xLx{0,1}[FYLIVMP]", "RxLE", options)

    @test result.matched
    @test result.query_relationship == "Degenerate Parent"
    @test result.search_relationship == "Variant Subsequence"
    @test result.matched_positions == 2
    @test result.match_ic ≈ 1.77 atol = 1e-2
    @test result.query_information ≈ 2.12 atol = 1e-2
    @test result.normalized_ic ≈ 0.835 atol = 1e-3
    @test result.score ≈ 1.669 atol = 1e-3
end
