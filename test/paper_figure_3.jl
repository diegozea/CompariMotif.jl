using TestItems

@testitem "paper Figure 3 examples" begin
    using Test
    using CompariMotif

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    # Edwards et al. (2008), Bioinformatics 24(10):1307-1309, Fig. 3.
    # Score literals below come from the figure's rounded display.
    # Figure 3 prints scores to two decimals; compare with a matching tolerance.
    cases = [
        (
            "R..S.P..L",
            "R[SFYW].S.P", # LIG_14-3-3_1
            "Degenerate Parent",
            "Variant Subsequence",
            2.54
        ),
        (
            "R..S.P..L",
            "[RHK][STALV].[ST].[PESRDIF]", # LIG_14-3-3_3
            "Complex Parent",
            "Complex Subsequence",
            2.37
        ),
        (
            "R..S.P..L",
            "[RK]..S[VI]..", # MOD_PK1
            "Complex Parent",
            "Complex Subsequence",
            1.40
        ),
        (
            "R..S.P..L",
            "[RK][RK].[ST]...", # MOD_PKA_1
            "Complex Parent",
            "Complex Subsequence",
            1.33
        ),
        (
            "FR..[ST]",
            "[RHK][STALV].[ST].[PESRDIF]", # LIG_14-3-3_3
            "Complex Overlap",
            "Complex Overlap",
            1.27
        ),
        (
            "FR..[ST]",
            "[RK][RK].[ST]...", # MOD_PKA_1
            "Complex Overlap",
            "Complex Overlap",
            1.33
        ),
        (
            "GR.[ST]..P",
            "[RHK][STALV].[ST].[PESRDIF]", # LIG_14-3-3_3
            "Complex Parent",
            "Complex Subsequence",
            0.89
        ),
        (
            "GR.[ST]..P",
            "R.[SYFWTQAD].[ST].[PLM]", # LIG_14-3-3_3
            "Complex Overlap",
            "Complex Overlap",
            0.97
        ),
        (
            "GR.[ST]..P",
            ".R.[ST]...", # MOD_PKA_2
            "Variant Match",
            "Degenerate Match",
            2.00
        )
    ]

    for (query, search, query_relationship, search_relationship, score) in cases
        @testset "$query vs $search" begin
            result = compare(query, search, options)
            @test result.matched
            @test result.query_relationship == query_relationship
            @test result.search_relationship == search_relationship
            @test result.score ≈ score atol = 0.01
        end
    end
end
