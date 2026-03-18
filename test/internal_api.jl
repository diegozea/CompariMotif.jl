@testitem "internal API docs pipeline example" begin
    # This test mirrors the worked example in `docs/src/internal_api.md`.
    # Keep this file and that page in sync when changing the example inputs,
    # the documented intermediate values, or the internal pipeline behavior.
    # The rounded numeric checks below intentionally mirror that docs walkthrough,
    # which is derived from the paper's Figure 1 example.

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    @testset "parse and normalize" begin
        parsed_query = CompariMotif._parse_motif("[KR]xLx{0,1}[FYLIVMP]", options)
        parsed_search = CompariMotif._parse_motif("RxLE", options)

        @test parsed_query.normalized == "[RK]xLx{0,1}[ILMFPYV]"
        @test parsed_search.normalized == "RxLE"
        @test length(parsed_query.alternatives) == 1
        @test length(parsed_search.alternatives) == 1

        spec = CompariMotif._alphabet_spec(options.alphabet)
        query_variants = CompariMotif._expand_variants(parsed_query, options, spec)
        search_variants = CompariMotif._expand_variants(parsed_search, options, spec)

        @test [variant.normalized for variant in query_variants] == [
            "[RK]xL[ILMFPYV]",
            "[RK]xLx[ILMFPYV]"
        ]
        @test round.([variant.information for variant in query_variants], digits = 2) ==
              [2.12, 2.12]
        @test only(search_variants).normalized == "RxLE"
        @test round(only(search_variants).information, digits = 3) == 3.0

        @testset "precise match pre-pass" begin
            found_precise,
            best_precise = CompariMotif._find_precise_match(
                query_variants,
                search_variants,
                options,
                spec
            )

            @test !found_precise
            @test best_precise === nothing
        end

        @testset "alignment scoring" begin
            candidate = CompariMotif._evaluate_alignment(
                query_variants[2],
                only(search_variants),
                0,
                options,
                spec
            )

            @test candidate.matched_pattern == "[rk]xLe"
            @test candidate.matched_positions == 2
            @test round(candidate.normalized_ic, digits = 3) == 0.835
        end
    end

    @testset "public result materialization" begin
        result = compare("[KR]xLx{0,1}[FYLIVMP]", "RxLE", options)

        @test result.matched
        @test result.query_relationship == "Degenerate Parent"
        @test result.search_relationship == "Variant Subsequence"
        @test round(result.match_ic, digits = 2) == 1.77
        @test round(result.score, digits = 3) == 1.669
    end
end
