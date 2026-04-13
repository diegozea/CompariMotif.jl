@testitem "docs regression guards" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    @testset "embedded grouped alternation follows oracle canonical order" begin
        grouped_overlap = compare("A(.|Q)L", ".[QL]", options)
        @test grouped_overlap.normalized_query == "(AQL)|(A.L)"
    end

    @testset "winner selection can be direction-dependent" begin
        forward = compare("A.{0,1}A.{2,3}\$", "A.KA", options)
        reverse = compare("A.KA", "A.{0,1}A.{2,3}\$", options)

        @test forward.query_relationship == "Complex Parent"
        @test forward.search_relationship == "Complex Subsequence"
        @test forward.matched_pattern == "Aaka"

        @test reverse.query_relationship == "Exact Overlap"
        @test reverse.search_relationship == "Exact Overlap"
        @test reverse.matched_pattern == "A"
    end

    @testset "contained-vs-overlap tie keeps only leading-prefix override" begin
        candidate = (; matched_positions,
            leading_exact_positions,
            overlap_length = 4,
            dual_wildcard_positions = 0,
            core_ic = 1.0,
            query_rel_length,
            search_rel_length) -> CompariMotif._Candidate(
            CompariMotif._MotifVariant(CompariMotif._Position[], "Q", 1.0),
            CompariMotif._MotifVariant(CompariMotif._Position[], "S", 1.0),
            CompariMotif._REL_COMPLEX,
            query_rel_length,
            CompariMotif._REL_COMPLEX,
            search_rel_length,
            "",
            matched_positions,
            1.0,
            1.0,
            core_ic,
            1.0,
            overlap_length,
            dual_wildcard_positions,
            0,
            leading_exact_positions
        )

        contained_same_prefix = candidate(;
            matched_positions = 2,
            leading_exact_positions = 1,
            overlap_length = 7,
            dual_wildcard_positions = 0,
            core_ic = 1.0,
            query_rel_length = CompariMotif._LEN_PARENT,
            search_rel_length = CompariMotif._LEN_SUBSEQUENCE
        )
        overlap_same_prefix = candidate(;
            matched_positions = 2,
            leading_exact_positions = 1,
            overlap_length = 4,
            dual_wildcard_positions = 2,
            core_ic = 1.2,
            query_rel_length = CompariMotif._LEN_OVERLAP,
            search_rel_length = CompariMotif._LEN_OVERLAP
        )
        @test !CompariMotif._is_better(contained_same_prefix, overlap_same_prefix)
        @test !CompariMotif._is_better(overlap_same_prefix, contained_same_prefix)

        contained_longer_prefix = candidate(;
            matched_positions = 2,
            leading_exact_positions = 2,
            overlap_length = 7,
            dual_wildcard_positions = 0,
            core_ic = 1.0,
            query_rel_length = CompariMotif._LEN_PARENT,
            search_rel_length = CompariMotif._LEN_SUBSEQUENCE
        )
        @test CompariMotif._is_better(contained_longer_prefix, overlap_same_prefix)
    end
end
