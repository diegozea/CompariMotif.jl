@testitem "internal motif normalization syntax" begin
    @test CompariMotif._normalize_motif("r[kR].{0,1}l") == "R[RK].{0,1}L"
    @test CompariMotif._normalize_motif("A[^P]x") == "A[ARNDCQEGHILKMFSTWYV]."
    @test CompariMotif._normalize_motif("[KKAQ]") == "[AQK]"
    @test_throws ArgumentError CompariMotif._normalize_motif("x(1,2)")
end

@testitem "alternation syntax" begin
    @test CompariMotif._normalize_motif("(rkli)") == "RKLI"
    @test CompariMotif._normalize_motif("(rkli)|(r[kr]l[iv])") == "(RKLI)|(R[RK]L[IV])"
    @test CompariMotif._normalize_motif("RKLI|R[KR]L[IV]") == "(RKLI)|(R[RK]L[IV])"
    @test CompariMotif._normalize_motif("A(K|Q)LI") == "(AKLI)|(AQLI)"
    @test CompariMotif._normalize_motif("(K|Q)") == "(K)|(Q)"
    @test CompariMotif._normalize_motif("(A|C){2}") == "(A{2})|(C{2})"
    @test CompariMotif._normalize_motif("(AC|GT){2}") == "(AC{2})|(GT{2})"
    @test CompariMotif._normalize_motif("R(KL)I") == "RKLI"
    @test CompariMotif._normalize_motif("(A.L)|(AQL)") == "(AQL)|(A.L)"
    @test CompariMotif._normalize_motif("A(.|Q)L") == "(AQL)|(A.L)"
    @test CompariMotif._normalize_motif("(Q[QL])|(Q.)") == "(Q.)|(Q[QL])"
    @test CompariMotif._normalize_motif("(.|Q)L") == "(QL)|(.L)"
    @test CompariMotif._normalize_motif("(^Q)|(Q)") == "(Q)|(^Q)"
    @test CompariMotif._normalize_motif("A((.|Q)|W)L") == "(AQL)|(AWL)|(A.L)"
    @test CompariMotif._normalize_motif("(.{0,1}AA)|(.AA)") == "(.AA)|(.{0,1}AA)"
    @test CompariMotif._normalize_motif("(QQ[ST])|(QQ{1,2})") == "(QQ{1,2})|(QQ[ST])"
    @test CompariMotif._normalize_motif("(A.{0,1})|(A{1,2})") == "(A{1,2})|(A.{0,1})"
    @test CompariMotif._normalize_motif("(A{1,4})|(A{1,3})") == "(A{1,3})|(A{1,4})"
    @test CompariMotif._normalize_motif("(A{1,3})|(A{2,3})") == "(A{1,3})|(A{2,3})"
    @test CompariMotif._normalize_motif("(A{2,3})|(A{1,3})") == "(A{1,3})|(A{2,3})"

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)
    # Oracle-derived canonical residue order for grouped alternation branches:
    # A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y.
    canonical_residues = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
        "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    residue_rank = Dict(residue => idx for (idx, residue) in enumerate(canonical_residues))

    expected_branch_order(a::String, b::String) = residue_rank[a] < residue_rank[b] ?
                                                  [a, b] : [b, a]

    result = compare("(RKLI)|(AQLI)", "AQLI", options)
    @test result.matched
    @test result.query_relationship == "Exact Match"
    @test result.normalized_query == "(AQLI)|(RKLI)"

    class_equivalence = compare("(K|Q)", "[KQ]", options)
    @test class_equivalence.matched
    @test class_equivalence.query_relationship == "Variant Match"
    @test class_equivalence.search_relationship == "Degenerate Match"

    redundant_grouping = compare("R(KL)I", "RKLI", options)
    @test redundant_grouping.matched
    @test redundant_grouping.query_relationship == "Exact Match"
    @test redundant_grouping.normalized_query == "RKLI"

    reordered = CompariMotif._expand_variants(
        CompariMotif._parse_motif("([ST]QQ)|(QQ[ST])", options),
        options,
        spec
    )
    @test getfield.(reordered, :normalized) == ["QQ[ST]", "[ST]QQ"]

    bounded_repeat = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(QQ[ST])|(QQ{1,2})", options),
        options,
        spec
    )
    @test getfield.(bounded_repeat, :normalized) == ["QQ", "QQQ", "QQ[ST]"]

    reverse_bounded_repeat = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(QQ{1,2})|(QQ[ST])", options),
        options,
        spec
    )
    @test getfield.(reverse_bounded_repeat, :normalized) == ["QQ", "QQQ", "QQ[ST]"]

    concrete_variant_tiebreak = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(A{1,2})|(A.{0,1})", options),
        options,
        spec
    )
    @test getfield.(concrete_variant_tiebreak, :normalized) == ["A", "AA", "A", "A."]

    reverse_concrete_variant_tiebreak = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(A.{0,1})|(A{1,2})", options),
        options,
        spec
    )
    @test getfield.(reverse_concrete_variant_tiebreak, :normalized) ==
          ["A", "AA", "A", "A."]

    overlapping_repeat_tiebreak = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(A{1,3})|(A{2,3})", options),
        options,
        spec
    )
    @test getfield.(overlapping_repeat_tiebreak, :normalized) ==
          ["A", "AA", "AAA", "AA", "AAA"]

    reverse_overlapping_repeat_tiebreak = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(A{2,3})|(A{1,3})", options),
        options,
        spec
    )
    @test getfield.(reverse_overlapping_repeat_tiebreak, :normalized) ==
          ["A", "AA", "AAA", "AA", "AAA"]

    for i in eachindex(canonical_residues), j in eachindex(canonical_residues)

        i == j && continue
        a = canonical_residues[i]
        b = canonical_residues[j]
        expected = expected_branch_order(a, b)

        plain = CompariMotif._expand_variants(
            CompariMotif._parse_motif("($a|$b)", options),
            options,
            spec
        )
        quantified = CompariMotif._expand_variants(
            CompariMotif._parse_motif("($a|$b){2}", options),
            options,
            spec
        )

        @test getfield.(plain, :normalized) == expected
        @test getfield.(quantified, :normalized) ==
              [repeat(expected[1], 2), repeat(expected[2], 2)]
    end

    mandatory_first = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(.AA)|(.{0,1}AA)", options),
        options,
        spec
    )
    @test getfield.(mandatory_first, :normalized) == [".AA", "AA", ".AA"]

    reverse_mandatory_first = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(.{0,1}AA)|(.AA)", options),
        options,
        spec
    )
    @test getfield.(reverse_mandatory_first, :normalized) == [".AA", "AA", ".AA"]

    repeat_tiebreak_probe = CompariMotif._expand_variants(
        CompariMotif._parse_motif("[DE].{0,2}[DE].{0,2}[ADE].{0,1}[CF]\$", options),
        options,
        spec
    )
    @test getfield.(repeat_tiebreak_probe, :normalized)[1:9] == [
        "[DE][DE][ADE][CF]\$",
        "[DE][DE][ADE].[CF]\$",
        "[DE][DE].[ADE][CF]\$",
        "[DE][DE].[ADE].[CF]\$",
        "[DE][DE]..[ADE][CF]\$",
        "[DE][DE]..[ADE].[CF]\$",
        "[DE].[DE][ADE][CF]\$",
        "[DE].[DE][ADE].[CF]\$",
        "[DE].[DE].[ADE][CF]\$"
    ]

    not_p = "[ARNDCQEGHILKMFSTWYV]"
    rk = "[RK]"
    multi_range_probe = CompariMotif._expand_variants(
        CompariMotif._parse_motif("[ST]{0,1}V.G[^P]{0,2}[KR]{1,2}", options),
        options,
        spec
    )
    @test getfield.(multi_range_probe, :normalized) == [
        "V.G" * rk,
        "V.G" * rk * rk,
        "V.G" * not_p * rk,
        "V.G" * not_p * rk * rk,
        "V.G" * not_p * not_p * rk,
        "V.G" * not_p * not_p * rk * rk,
        "[ST]V.G" * rk,
        "[ST]V.G" * rk * rk,
        "[ST]V.G" * not_p * rk,
        "[ST]V.G" * not_p * rk * rk,
        "[ST]V.G" * not_p * not_p * rk,
        "[ST]V.G" * not_p * not_p * rk * rk
    ]

    bounded_repeat_match = compare("(QQ[ST])|(QQ{1,2})", "QQ", options)
    @test bounded_repeat_match.matched
    @test bounded_repeat_match.query_relationship == "Exact Match"
    @test bounded_repeat_match.search_relationship == "Exact Match"
    @test bounded_repeat_match.matched_pattern == "QQ"

    mandatory_repeat_match = compare("(.AA)|(.{0,1}AA)", "AA", options)
    @test mandatory_repeat_match.matched
    @test mandatory_repeat_match.query_relationship == "Exact Parent"
    @test mandatory_repeat_match.search_relationship == "Exact Subsequence"
    @test mandatory_repeat_match.matched_pattern == "AA"

    explicit_tie_branch_order = compare("(A.QQ)|(QQA)", "(.QQ)|(.{0,1}QQ)", options)
    @test explicit_tie_branch_order.matched
    @test explicit_tie_branch_order.query_relationship == "Exact Parent"
    @test explicit_tie_branch_order.search_relationship == "Exact Subsequence"
    @test explicit_tie_branch_order.matched_pattern == ".QQ"

    overlapping_repeat_match = compare("(A{1,3}|A{2,3})", "[QA]", options)
    @test overlapping_repeat_match.matched
    @test overlapping_repeat_match.query_relationship == "Variant Match"
    @test overlapping_repeat_match.search_relationship == "Degenerate Match"
    @test overlapping_repeat_match.matched_pattern == "[aq]"

    embedded_group_variants = CompariMotif._expand_variants(
        CompariMotif._parse_motif("A(Q|.)L", options),
        options,
        spec
    )
    @test getfield.(embedded_group_variants, :normalized) == ["AQL", "A.L"]

    reverse_embedded_group_variants = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(A.L)|(AQL)", options),
        options,
        spec
    )
    @test getfield.(reverse_embedded_group_variants, :normalized) == ["AQL", "A.L"]

    embedded_specificity_variants = CompariMotif._expand_variants(
        CompariMotif._parse_motif("Q(.|[QL])", options),
        options,
        spec
    )
    @test getfield.(embedded_specificity_variants, :normalized) == ["Q.", "Q[QL]"]

    reverse_embedded_specificity_variants = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(Q[QL])|(Q.)", options),
        options,
        spec
    )
    @test getfield.(reverse_embedded_specificity_variants, :normalized) == ["Q.", "Q[QL]"]

    anchor_branch_variants = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(^Q)|(Q)", options),
        options,
        spec
    )
    @test getfield.(anchor_branch_variants, :normalized) == ["Q", "^Q"]

    embedded_specificity_match = compare("Q(.|[QL])", "[QL].", options)
    @test embedded_specificity_match.matched
    @test embedded_specificity_match.query_relationship == "Variant Match"
    @test embedded_specificity_match.search_relationship == "Degenerate Match"

    embedded_specificity_overlap = compare("Q(.|[QL])", ".[QL]", options)
    @test embedded_specificity_overlap.matched
    @test embedded_specificity_overlap.query_relationship == "Variant Overlap"
    @test embedded_specificity_overlap.search_relationship == "Degenerate Overlap"

    nested_group_variants = CompariMotif._expand_variants(
        CompariMotif._parse_motif("A((.|Q)|W)L", options),
        options,
        spec
    )
    @test getfield.(nested_group_variants, :normalized) == ["AQL", "AWL", "A.L"]
end

@testitem "_branch_order_tokens comment examples" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)

    branch_and_order(motif::AbstractString) = begin
        tokens = CompariMotif._parse_linear_tokens(motif, spec)
        return tokens, CompariMotif._branch_order_tokens(tokens)
    end

    tokens_r2, order_r2 = branch_and_order("R{2}")
    @test [(tok.min_repeat, tok.max_repeat) for tok in order_r2] == [(1, 1), (1, 1)]
    @test all(tok -> CompariMotif._same_position(tok.position, tokens_r2[1].position), order_r2)

    tokens_r24, order_r24 = branch_and_order("R{2,4}")
    @test [(tok.min_repeat, tok.max_repeat) for tok in order_r24] ==
          [(1, 1), (1, 1), (0, 2)]
    @test all(tok -> CompariMotif._same_position(tok.position, tokens_r24[1].position), order_r24)

    tokens_r03, order_r03 = branch_and_order("R{0,3}")
    @test [(tok.min_repeat, tok.max_repeat) for tok in order_r03] == [(0, 3)]
    @test all(tok -> CompariMotif._same_position(tok.position, tokens_r03[1].position), order_r03)
end

@testitem "_order_token_repeat_key docstring examples" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)

    exact = only(CompariMotif._parse_linear_tokens("R{3}", spec))
    ranged = only(CompariMotif._parse_linear_tokens("R{2,4}", spec))
    optional_tail = last(CompariMotif._branch_order_tokens([ranged]))

    @test CompariMotif._order_token_repeat_key(exact) == (0, 3)
    @test CompariMotif._order_token_repeat_key(ranged) == (1, 4)
    @test CompariMotif._order_token_repeat_key(optional_tail) == (1, 2)
end

@testitem "branch ordering fallback helpers" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)

    left_repeat_overlap = CompariMotif._branch_order_tokens(
        CompariMotif._parse_linear_tokens("CA{0,1}", spec)
    )
    right_repeat_overlap = CompariMotif._branch_order_tokens(
        CompariMotif._parse_linear_tokens("AA", spec)
    )
    @test CompariMotif._repeat_overlap_ordering(
        left_repeat_overlap,
        right_repeat_overlap,
        2
    ) === nothing

    tied_tokens = CompariMotif._parse_linear_tokens("AA", spec)
    tied_order_tokens = CompariMotif._branch_order_tokens(tied_tokens)
    earlier_branch = (;
        branch_index = 1,
        tokens = tied_tokens,
        order_tokens = tied_order_tokens
    )
    later_branch = (;
        branch_index = 2,
        tokens = copy(tied_tokens),
        order_tokens = copy(tied_order_tokens)
    )
    @test CompariMotif._branch_tokens_isless(earlier_branch, later_branch, spec)
    @test !CompariMotif._branch_tokens_isless(later_branch, earlier_branch, spec)
end

@testitem "oracle parser corner cases" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)

    @test CompariMotif._normalize_motif(" AC ") == "AC"
    @test CompariMotif._normalize_motif("A C") == "A"
    @test CompariMotif._normalize_motif("A (K|Q) L I") == "A"
    @test CompariMotif._normalize_motif("A)") == "A"
    @test CompariMotif._normalize_motif("A{-1}") == "A"
    @test CompariMotif._normalize_motif("A{0}") == "A"
    @test CompariMotif._normalize_motif("A{0}CD") == "ACD"
    @test CompariMotif._normalize_motif("^{2}A") == "^{2}A"
    @test CompariMotif._normalize_motif("A\${2}") == "A\${2}"
    @test CompariMotif._normalize_motif("A{0,1}") == "A{0,1}"
    @test CompariMotif._normalize_motif("A{2,1}") == "A{2,1}"
    @test CompariMotif._normalize_motif("(Q|.)") == "(Q)|(.)"
    @test CompariMotif._normalize_motif("(Q|x)") == "(Q)|(.)"
    @test CompariMotif._normalize_motif("(Q|X)") == "(Q)|(.)"
    @test CompariMotif._normalize_motif("(Q|[ARNDCQEGHILKMFPSTWYV])") == "(Q)|(.)"

    retained = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("A{0}", options), options, spec))
    @test retained.normalized == "A"

    grouped = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(A|C){2}", options), options, spec)
    @test sort(getfield.(grouped, :normalized)) == ["AA", "CC"]
    grouped_long = CompariMotif._expand_variants(
        CompariMotif._parse_motif("(AC|GT){2}", options), options, spec)
    @test sort(getfield.(grouped_long, :normalized)) == ["ACC", "GTT"]

    anchor_repeat = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("^{2}A", options), options, spec))
    @test anchor_repeat.normalized == "^^A"

    impossible = CompariMotif._expand_variants(
        CompariMotif._parse_motif("A{2,1}", options), options, spec)
    @test isempty(impossible)
    @test !compare("A{2,1}", "A", options).matched

    dot_group = CompariMotif._parse_motif("(Q|.)", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in dot_group.alternatives] == ["Q", "."]
    @test compare("(Q|.)", "Q", options).matched
    reverse_dot_group = CompariMotif._parse_motif("(.|Q)", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in reverse_dot_group.alternatives] == ["Q", "."]
    @test compare("(.|Q)", "Q", options).matched

    upper_x_group = CompariMotif._parse_motif("(Q|X)", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in upper_x_group.alternatives] == ["Q", "."]
    @test compare("(Q|X)", "Q", options).matched
    @test compare("Q|X", "Q", options).matched

    lower_x_group = CompariMotif._parse_motif("(Q|x)", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in lower_x_group.alternatives] == ["Q", "."]
    @test compare("(Q|x)", "Q", options).matched

    quantified_wildcard_group = CompariMotif._parse_motif("(Q|.){2}", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in quantified_wildcard_group.alternatives] == ["Q{2}", ".{2}"]
    @test compare("(Q|.){2}", "QQ", options).matched
    @test compare("(Q|x){2}", "QQ", options).matched
    @test compare("(Q|X){2}", "QQ", options).matched

    embedded_dot_group = CompariMotif._parse_motif("A(Q|.)L", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in embedded_dot_group.alternatives] ==
          ["AQL", "A.L"]
    @test compare("A(Q|.)L", "AQL", options).matched

    embedded_specificity_group = CompariMotif._parse_motif("Q(.|[QL])", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in embedded_specificity_group.alternatives] ==
          ["Q.", "Q[QL]"]
    reverse_embedded_specificity_group = CompariMotif._parse_motif("(Q[QL])|(Q.)", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in reverse_embedded_specificity_group.alternatives] ==
          ["Q.", "Q[QL]"]

    anchor_group = CompariMotif._parse_motif("(^Q)|(Q)", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in anchor_group.alternatives] ==
          ["Q", "^Q"]

    suffix_dot_group = CompariMotif._parse_motif("(Q|.)L", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in suffix_dot_group.alternatives] ==
          ["QL", ".L"]
    @test compare("(Q|.)L", "QL", options).matched

    preserved_suffix_group = CompariMotif._parse_motif("(.|Q)L", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in preserved_suffix_group.alternatives] ==
          ["QL", ".L"]
    preserved_suffix_quantifier_group = CompariMotif._parse_motif("(.|Q)C{2}", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in preserved_suffix_quantifier_group.alternatives] ==
          ["QC{2}", ".C{2}"]

    prefix_dot_group = CompariMotif._parse_motif("A(Q|.)", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in prefix_dot_group.alternatives] ==
          ["AQ", "A."]
    @test compare("A(Q|.)", "AQ", options).matched

    nested_group = CompariMotif._parse_motif("A((.|Q)|W)L", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in nested_group.alternatives] ==
          ["AQL", "AWL", "A.L"]

    embedded_upper_x_group = CompariMotif._parse_motif("A(Q|X)L", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in embedded_upper_x_group.alternatives] ==
          ["AQL", "A.L"]
    @test compare("A(Q|X)L", "AQL", options).matched

    embedded_lower_x_group = CompariMotif._parse_motif("A(Q|x)L", options)
    @test [CompariMotif._normalized_from_tokens(tokens)
           for tokens in embedded_lower_x_group.alternatives] ==
          ["AQL", "A.L"]
    @test compare("A(Q|x)L", "AQL", options).matched

    retained_full_alphabet_group = compare("(Q|[ARNDCQEGHILKMFPSTWYV])", "Q", options)
    @test retained_full_alphabet_group.matched
    @test retained_full_alphabet_group.normalized_query == "(Q)|(.)"
    @test compare(retained_full_alphabet_group.normalized_query, "Q", options).matched

    dot_result = compare("(Q|.)", "Q", options)
    lower_x_result = compare("(Q|x)", "Q", options)
    upper_x_result = compare("(Q|X)", "Q", options)
    for result in (dot_result, lower_x_result, upper_x_result)
        @test result.matched
        @test result.query_relationship == "Exact Match"
        @test result.search_relationship == "Exact Match"
        @test result.matched_pattern == "Q"
        @test result.normalized_query == "(Q)|(.)"
    end

    @test_throws ArgumentError CompariMotif._normalize_motif("A[Q")
    @test_throws ArgumentError CompariMotif._normalize_motif("A{-1,2}")
end

@testitem "wildcard token equivalence" begin
    @test CompariMotif._normalize_motif("A.Xx"; alphabet = ProteinAlphabet()) == "A..."
    @test CompariMotif._normalize_motif("A.Xx"; alphabet = DNAAlphabet()) == "A..."
    @test CompariMotif._normalize_motif("A.Xx"; alphabet = RNAAlphabet()) == "A..."
end

@testitem "position information content" begin
    for alphabet in (ProteinAlphabet(), DNAAlphabet(), RNAAlphabet())
        options = ComparisonOptions(; alphabet)
        spec = CompariMotif._alphabet_spec(alphabet)

        wildcard = only(CompariMotif._expand_variants(
            CompariMotif._parse_motif("x", options), options, spec))
        fixed = only(CompariMotif._expand_variants(
            CompariMotif._parse_motif("A", options), options, spec))

        @test wildcard.information == 0.0
        @test CompariMotif._position_ic(only(wildcard.positions), spec) == 0.0

        @test fixed.information ≈ 1.0 atol = 1e-12
        @test CompariMotif._position_ic(only(fixed.positions), spec) ≈ 1.0 atol = 1e-12
    end
end

@testitem "positive character classes ignore duplicates" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    duplicate_class = compare("[KKAQ]", "Q", options)
    unique_class = compare("[AQK]", "Q", options)
    double_a = compare("[AA]", "A", options)
    triple_a = compare("[AAA]", "A", options)
    branch_selection = compare("N{1,2}Q([KKAQ]|[NW]{1,2}[TI])", "NQ[QK]", options)

    @test CompariMotif._normalize_motif("[AA]") == "A"
    @test CompariMotif._normalize_motif("[AAA]") == "A"
    @test CompariMotif._normalize_motif("[KKAQ]") == "[AQK]"

    @test duplicate_class.match_ic ≈ unique_class.match_ic atol = 1e-12
    @test duplicate_class.query_information ≈ unique_class.query_information atol = 1e-12
    @test duplicate_class.query_relationship == unique_class.query_relationship
    @test duplicate_class.search_relationship == unique_class.search_relationship
    @test duplicate_class.matched_pattern == unique_class.matched_pattern

    @test double_a.query_relationship == "Exact Match"
    @test double_a.search_relationship == "Exact Match"
    @test double_a.match_ic ≈ 1.0 atol = 1e-12
    @test triple_a.query_relationship == "Exact Match"
    @test triple_a.search_relationship == "Exact Match"
    @test triple_a.match_ic ≈ 1.0 atol = 1e-12

    @test branch_selection.query_relationship == "Degenerate Match"
    @test branch_selection.search_relationship == "Variant Match"
    @test branch_selection.match_ic ≈ 2.6332742086579152 atol = 1e-12
    @test branch_selection.normalized_ic ≈ 1.0 atol = 1e-12
    @test branch_selection.core_ic ≈ 0.9511137350628183 atol = 1e-12
    @test branch_selection.score ≈ 3.0 atol = 1e-12
    @test branch_selection.query_information ≈ 2.6332742086579152 atol = 1e-12
    @test branch_selection.normalized_query == "(N{1,2}Q[AQK])|(N{1,2}Q[NW]{1,2}[IT])"
end

@testitem "position comparison primitives" begin
    using CompariMotif: CompariMotif, ComparisonOptions

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0, mismatches = 1)
    spec = CompariMotif._alphabet_spec(options.alphabet)
    residue_frequencies = CompariMotif._residue_frequency_vector(options, spec)

    function first_position(pattern::AbstractString)
        variant = only(CompariMotif._expand_variants(
            CompariMotif._parse_motif(pattern, options), options, spec))
        return only(variant.positions)
    end

    wildcard_exact = CompariMotif._compare_positions(
        first_position("x"),
        first_position("x"),
        options,
        spec,
        residue_frequencies
    )
    @test !wildcard_exact.hard_mismatch
    @test !wildcard_exact.mismatch
    @test wildcard_exact.relation == CompariMotif._REL_EXACT
    @test wildcard_exact.ic == 0.0
    @test wildcard_exact.core_ic_denominator == 0.0
    @test !wildcard_exact.contributes_position

    wildcard_degenerate = CompariMotif._compare_positions(
        first_position("x"),
        first_position("A"),
        options,
        spec,
        residue_frequencies
    )
    @test !wildcard_degenerate.hard_mismatch
    @test !wildcard_degenerate.mismatch
    @test wildcard_degenerate.relation == CompariMotif._REL_DEGENERATE
    @test wildcard_degenerate.ic == 0.0
    @test wildcard_degenerate.core_ic_denominator ≈ 1.0 atol = 1e-12
    @test !wildcard_degenerate.contributes_position

    wildcard_variant = CompariMotif._compare_positions(
        first_position("A"),
        first_position("x"),
        options,
        spec,
        residue_frequencies
    )
    @test !wildcard_variant.hard_mismatch
    @test !wildcard_variant.mismatch
    @test wildcard_variant.relation == CompariMotif._REL_VARIANT
    @test wildcard_variant.ic == 0.0
    @test wildcard_variant.core_ic_denominator ≈ 1.0 atol = 1e-12
    @test !wildcard_variant.contributes_position

    nterm_exact = CompariMotif._compare_positions(
        first_position(raw"^"),
        first_position(raw"^"),
        options,
        spec,
        residue_frequencies
    )
    @test !nterm_exact.hard_mismatch
    @test !nterm_exact.mismatch
    @test nterm_exact.relation == CompariMotif._REL_EXACT
    @test nterm_exact.ic == 1.0
    @test nterm_exact.core_ic_denominator == 1.0
    @test nterm_exact.contributes_position

    cterm_exact = CompariMotif._compare_positions(
        first_position(raw"$"),
        first_position(raw"$"),
        options,
        spec,
        residue_frequencies
    )
    @test !cterm_exact.hard_mismatch
    @test !cterm_exact.mismatch
    @test cterm_exact.relation == CompariMotif._REL_EXACT
    @test cterm_exact.ic == 1.0
    @test cterm_exact.core_ic_denominator == 1.0
    @test cterm_exact.contributes_position

    opposite_anchor = CompariMotif._compare_positions(
        first_position(raw"^"),
        first_position(raw"$"),
        options,
        spec,
        residue_frequencies
    )
    @test opposite_anchor.hard_mismatch
    @test opposite_anchor.mismatch
    @test opposite_anchor.relation == CompariMotif._REL_COMPLEX
    @test opposite_anchor.ic == 0.0
    @test opposite_anchor.core_ic_denominator == 0.0
    @test !opposite_anchor.contributes_position
end

@testitem "alignment helper edge coverage" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)

    function only_variant(pattern::AbstractString)
        parsed = CompariMotif._parse_motif(pattern, options)
        return only(CompariMotif._expand_variants(parsed, options, spec))
    end

    anchor_error = try
        CompariMotif._anchor_symbol(
            CompariMotif._Position(CompariMotif._RESIDUE, zero(CompariMotif.ResidueMask))
        )
        nothing
    catch err
        err
    end
    @test anchor_error isa ErrorException
    if anchor_error isa ErrorException
        @test occursin("Expected anchor position", sprint(showerror, anchor_error))
    end

    exact_query = only_variant("AK")
    exact_search = only_variant("AK")
    @test CompariMotif._is_exact_contained_alignment(
        exact_query.positions,
        exact_search.positions,
        1,
        1,
        length(exact_query.positions),
        false,
        false,
        false,
        false,
        spec.mask
    )

    wildcard_only = only_variant("...")
    @test CompariMotif._first_informative_index(wildcard_only.positions, spec.mask) ===
          nothing
end

@testitem "wildcard positions do not count toward matched positions" begin
    using CompariMotif: ComparisonOptions, compare

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    leading_wildcard = compare("AK", "A.", options)
    @test leading_wildcard.matched
    @test leading_wildcard.matched_positions == 1
    @test leading_wildcard.match_ic ≈ 1.0 atol = 1e-12
    @test leading_wildcard.core_ic ≈ 0.5 atol = 1e-12

    trailing_wildcard = compare("A.", "AK", options)
    @test trailing_wildcard.matched
    @test trailing_wildcard.matched_positions == 1
    @test trailing_wildcard.match_ic ≈ 1.0 atol = 1e-12
    @test trailing_wildcard.core_ic ≈ 0.5 atol = 1e-12

    wildcard_only = compare("x", "A", options)
    @test !wildcard_only.matched
end

@testitem "non-uniform residue frequencies" begin
    cases = (
        (
            ProteinAlphabet(),
            Dict(aa => 1.0 for aa in "ARNDCQEGHILKMFPSTWYV"),
            'A',
            'C',
            "[AC]"
        ),
        (
            DNAAlphabet(),
            Dict('A' => 7.0, 'C' => 2.0, 'G' => 1.0, 'T' => 1.0),
            'A',
            'C',
            "[AC]"
        ),
        (
            RNAAlphabet(),
            Dict('A' => 7.0, 'C' => 2.0, 'G' => 1.0, 'U' => 1.0),
            'A',
            'C',
            "[AC]"
        )
    )

    for (alphabet, weights, common, rare, class_pattern) in cases
        if alphabet isa ProteinAlphabet
            weights[common] = 70.0
            weights[rare] = 20.0
        end
        options = ComparisonOptions(;
            alphabet,
            residue_frequencies = weights,
            min_shared_positions = 1,
            normalized_ic_cutoff = 0.0
        )
        spec = CompariMotif._alphabet_spec(alphabet)
        vector = CompariMotif._residue_frequency_vector(options, spec)

        common_variant = only(CompariMotif._expand_variants(
            CompariMotif._parse_motif(string(common), options), options, spec))
        class_variant = only(CompariMotif._expand_variants(
            CompariMotif._parse_motif(class_pattern, options), options, spec))

        common_expected = -log(options.residue_frequencies[common]) /
                          log(length(spec.chars))
        class_expected = -log(
            options.residue_frequencies['A'] + options.residue_frequencies['C']
        ) / log(length(spec.chars))

        @test CompariMotif._position_ic(only(common_variant.positions), spec, vector) ≈
              common_expected atol = 1e-12
        @test common_variant.information ≈ common_expected atol = 1e-12
        @test class_variant.information ≈ class_expected atol = 1e-12

        result = compare(string(common), class_pattern, options)
        @test result.matched
        @test result.match_ic ≈ class_expected atol = 1e-12
        @test result.query_information ≈ common_expected atol = 1e-12
        @test result.search_information ≈ class_expected atol = 1e-12
    end

    protein_weights = Dict(
        'A' => 0.7, 'C' => 0.2, 'D' => 0.01, 'E' => 0.01, 'F' => 0.005,
        'G' => 0.005, 'H' => 0.005, 'I' => 0.005, 'K' => 0.005, 'L' => 0.005,
        'M' => 0.005, 'N' => 0.005, 'P' => 0.005, 'Q' => 0.005, 'R' => 0.005,
        'S' => 0.005, 'T' => 0.005, 'V' => 0.005, 'W' => 0.005, 'Y' => 0.005
    )
    weighted_options = ComparisonOptions(;
        residue_frequencies = protein_weights,
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0
    )
    weighted_spec = CompariMotif._alphabet_spec(weighted_options.alphabet)
    weighted_union = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("[ACD]", weighted_options), weighted_options, weighted_spec))
    weighted_query = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("[AD]", weighted_options), weighted_options, weighted_spec))
    weighted_search = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("[AC]", weighted_options), weighted_options, weighted_spec))

    weighted_complex = compare("[AD]", "[AC]", weighted_options)
    @test weighted_complex.matched
    @test weighted_complex.query_relationship == "Complex Match"
    @test weighted_complex.search_relationship == "Complex Match"
    @test weighted_complex.match_ic ≈ weighted_union.information atol = 1e-12
    @test weighted_complex.normalized_ic ≈
          weighted_union.information /
          min(weighted_query.information, weighted_search.information) atol = 1e-12
end

@testitem "precise match phase precedes sliding window" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)

    query = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("RKLI", options), options, spec))
    search = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("ARKLIT", options), options, spec))

    found_precise,
    candidate = CompariMotif._find_precise_match([query], [search], options, spec)

    @test found_precise
    @test candidate !== nothing
    @test candidate.query_relationship_type == CompariMotif._REL_EXACT
    @test candidate.query_relationship_length == CompariMotif._LEN_SUBSEQUENCE
    @test candidate.search_relationship_type == CompariMotif._REL_EXACT
    @test candidate.search_relationship_length == CompariMotif._LEN_PARENT

    # Exact wildcard equality is detected in the precise-match phase, but can
    # still fail downstream thresholds because wildcards contribute no IC.
    strict_options = ComparisonOptions()
    strict_spec = CompariMotif._alphabet_spec(strict_options.alphabet)
    wildcard = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("x", strict_options), strict_options, strict_spec))

    found_wildcard_precise,
    wildcard_candidate = CompariMotif._find_precise_match(
        [wildcard], [wildcard], strict_options, strict_spec)

    @test found_wildcard_precise
    @test wildcard_candidate === nothing

    alternation_result = compare("x|RKL", "xKL", ComparisonOptions())
    @test alternation_result.matched
    @test alternation_result.query_relationship == "Variant Match"
    @test alternation_result.search_relationship == "Degenerate Match"

    # Precise exact-subsequence hits from one expanded branch should seed the
    # search order, but not suppress stronger candidates from another branch.
    branch_competition = compare("R|RKLI", "R[KR]L[IV]", options)
    @test branch_competition.matched
    @test branch_competition.query_relationship == "Variant Match"
    @test branch_competition.search_relationship == "Degenerate Match"
    @test branch_competition.matched_pattern == "R[rk]L[iv]"
    @test branch_competition.matched_positions == 4

    # When match_ic and matched_positions tie across expanded branches, the
    # higher score should win even if it was encountered later in search order.
    # The numeric winners below intentionally mirror the committed oracle score-
    # tiebreak fixture (`data/fixtures/oracle_score_tiebreak_probe_normalized.tsv`),
    # but stay inline because this unit test is meant to exercise the internal
    # tie-break path without depending on fixture I/O.
    zero_fixed_tiebreak = compare("A", "[AC][AC]|[AC]", options)
    @test zero_fixed_tiebreak.matched
    @test zero_fixed_tiebreak.query_relationship == "Variant Match"
    @test zero_fixed_tiebreak.search_relationship == "Degenerate Match"
    @test zero_fixed_tiebreak.score ≈ 1.0 atol = 1e-12
    @test zero_fixed_tiebreak.search_information ≈ 0.7686217868402407 atol = 1e-12

    nonzero_fixed_tiebreak = compare("AC", "AA|[AC]A", options)
    @test nonzero_fixed_tiebreak.matched
    @test nonzero_fixed_tiebreak.query_relationship == "Exact Overlap"
    @test nonzero_fixed_tiebreak.search_relationship == "Exact Overlap"
    @test nonzero_fixed_tiebreak.score ≈ 0.5654120103239064 atol = 1e-12
    @test nonzero_fixed_tiebreak.search_information ≈ 1.7686217868402407 atol = 1e-12

    shifted_nonzero_fixed_tiebreak = compare("CA", "AA|A[AC]", options)
    @test shifted_nonzero_fixed_tiebreak.matched
    @test shifted_nonzero_fixed_tiebreak.query_relationship == "Exact Overlap"
    @test shifted_nonzero_fixed_tiebreak.search_relationship == "Exact Overlap"
    @test shifted_nonzero_fixed_tiebreak.score ≈ 0.5654120103239064 atol = 1e-12
    @test shifted_nonzero_fixed_tiebreak.search_information ≈ 1.7686217868402407 atol = 1e-12

    subsequence_tiebreak = compare("AC", "AA|A", options)
    @test subsequence_tiebreak.matched
    @test subsequence_tiebreak.query_relationship == "Exact Parent"
    @test subsequence_tiebreak.search_relationship == "Exact Subsequence"
    @test subsequence_tiebreak.score ≈ 1.0 atol = 1e-12
    @test subsequence_tiebreak.search_information ≈ 1.0 atol = 1e-12

    oracle_shift_tiebreak = compare("AA", "Ax", options)
    @test oracle_shift_tiebreak.matched
    @test oracle_shift_tiebreak.query_relationship == "Exact Overlap"
    @test oracle_shift_tiebreak.search_relationship == "Exact Overlap"
    @test oracle_shift_tiebreak.matched_pattern == "A"

    reverse_oracle_shift_tiebreak = compare("Ax", "AA", options)
    @test reverse_oracle_shift_tiebreak.matched
    @test reverse_oracle_shift_tiebreak.query_relationship == "Degenerate Match"
    @test reverse_oracle_shift_tiebreak.search_relationship == "Variant Match"
    @test reverse_oracle_shift_tiebreak.matched_pattern == "Aa"

    less_ambiguous_tiebreak = compare("DE.", "[DE]{0,2}[DE].{0,2}[DE].{2,4}[CF].[CF]", options)
    @test less_ambiguous_tiebreak.matched
    @test less_ambiguous_tiebreak.query_relationship == "Variant Subsequence"
    @test less_ambiguous_tiebreak.search_relationship == "Degenerate Parent"
    @test less_ambiguous_tiebreak.matched_pattern == "[de][de]."

    fuller_coverage_tiebreak = compare(".(A)..C\$", "[DE].{0,2}[DE].{0,2}[ADE].{0,1}[CF]\$", options)
    @test fuller_coverage_tiebreak.matched
    @test fuller_coverage_tiebreak.query_relationship == "Complex Match"
    @test fuller_coverage_tiebreak.search_relationship == "Complex Match"
    @test fuller_coverage_tiebreak.matched_pattern == "[de]a[de][ade][cf]\$"

    shorter_complex_overlap_tiebreak = compare("A.{0,1}A.{2,3}\$", "[AS].{0,2}S\$", options)
    @test shorter_complex_overlap_tiebreak.matched
    @test shorter_complex_overlap_tiebreak.query_relationship == "Complex Parent"
    @test shorter_complex_overlap_tiebreak.search_relationship == "Complex Subsequence"
    @test shorter_complex_overlap_tiebreak.matched_pattern == "[as].s\$"

    single_position_overlap_wins = compare("AA.", "[AS].{0,2}S\$", options)
    @test single_position_overlap_wins.matched
    @test single_position_overlap_wins.query_relationship == "Complex Overlap"
    @test single_position_overlap_wins.search_relationship == "Complex Overlap"
    @test single_position_overlap_wins.matched_pattern == "[as]s"

    contained_complex_tiebreak = compare("[ST]{0,1}V.G[^P]{0,2}[KR]{1,2}", "R...[KR]R.", options)
    @test contained_complex_tiebreak.matched
    @test contained_complex_tiebreak.query_relationship == "Complex Subsequence"
    @test contained_complex_tiebreak.search_relationship == "Complex Parent"
    @test contained_complex_tiebreak.matched_positions == 2
    @test contained_complex_tiebreak.score ≈ 1.1104756749277043 atol = 1e-12

    single_position_subsequence_falls_back_to_overlap = compare(
        "R.x{0,1}L",
        "H[LM]H(([KR][^H].)|(.[^H][KR]))",
        options
    )
    @test single_position_subsequence_falls_back_to_overlap.matched
    @test single_position_subsequence_falls_back_to_overlap.query_relationship ==
          "Complex Overlap"
    @test single_position_subsequence_falls_back_to_overlap.search_relationship ==
          "Complex Overlap"
    @test single_position_subsequence_falls_back_to_overlap.matched_pattern == "h[lm]"

    contained_complex_overlap_regression = compare(
        "[KR]{1,4}[KR].[KR]W.",
        "[^DE]((K[RK])|(RK))[KRP][KR][^DE]",
        options
    )
    @test contained_complex_overlap_regression.matched
    @test contained_complex_overlap_regression.query_relationship ==
          "Complex Parent"
    @test contained_complex_overlap_regression.search_relationship ==
          "Complex Subsequence"
    @test contained_complex_overlap_regression.matched_pattern ==
          "[arncqghilkmfpstwyv][rk][RK][rkp][RK][arncqghilkmfpstwyv]"

    single_position_subsequence_can_still_beat_overlap = compare(
        "R.x{0,1}L",
        "[KR]{1,4}[KR].[KR]W.",
        options
    )
    @test single_position_subsequence_can_still_beat_overlap.matched
    @test single_position_subsequence_can_still_beat_overlap.query_relationship ==
          "Complex Subsequence"
    @test single_position_subsequence_can_still_beat_overlap.search_relationship ==
          "Complex Parent"
    @test single_position_subsequence_can_still_beat_overlap.matched_pattern == "[rk][rk]l"

    reverse_single_position_subsequence_can_still_beat_overlap = compare(
        "[KR]{1,4}[KR].[KR]W.",
        "R.x{0,1}L",
        options
    )
    @test reverse_single_position_subsequence_can_still_beat_overlap.matched
    @test reverse_single_position_subsequence_can_still_beat_overlap.query_relationship ==
          "Complex Parent"
    @test reverse_single_position_subsequence_can_still_beat_overlap.search_relationship ==
          "Complex Subsequence"
    @test reverse_single_position_subsequence_can_still_beat_overlap.matched_pattern ==
          "[rk]wl"

    repeat_expansion_overlap_tie_preserves_first_variant = compare(
        "[KR]{1,4}[KR].[KR]W.",
        "[NTSRM]..[TMSRG]Q.R..",
        options
    )
    @test repeat_expansion_overlap_tie_preserves_first_variant.matched
    @test repeat_expansion_overlap_tie_preserves_first_variant.query_relationship ==
          "Complex Overlap"
    @test repeat_expansion_overlap_tie_preserves_first_variant.search_relationship ==
          "Complex Overlap"
    @test repeat_expansion_overlap_tie_preserves_first_variant.matched_pattern ==
          "[rk][rk].[rk]"

    reverse_repeat_expansion_overlap_tie_preserves_first_variant = compare(
        "[NTSRM]..[TMSRG]Q.R..",
        "[KR]{1,4}[KR].[KR]W.",
        options
    )
    @test reverse_repeat_expansion_overlap_tie_preserves_first_variant.matched
    @test reverse_repeat_expansion_overlap_tie_preserves_first_variant.query_relationship ==
          "Complex Overlap"
    @test reverse_repeat_expansion_overlap_tie_preserves_first_variant.search_relationship ==
          "Complex Overlap"
    @test reverse_repeat_expansion_overlap_tie_preserves_first_variant.matched_pattern ==
          "[rk][rk]."

    wildcard_padding_subsequence_falls_back_to_overlap = compare(
        "[^DE]((K[RK])|(RK))[KRP][KR][^DE]",
        "[KR]{1,4}[KR].[KR]W.",
        options
    )
    @test wildcard_padding_subsequence_falls_back_to_overlap.matched
    @test wildcard_padding_subsequence_falls_back_to_overlap.query_relationship ==
          "Complex Overlap"
    @test wildcard_padding_subsequence_falls_back_to_overlap.search_relationship ==
          "Complex Overlap"
    @test wildcard_padding_subsequence_falls_back_to_overlap.matched_pattern ==
          "[rk][RK][rkp][RK][arncqghilkmfpstwyv]"
    @test wildcard_padding_subsequence_falls_back_to_overlap.core_ic ≈
          0.7299838398707831 atol = 1e-12

    contained_complex_review_regression = compare(
        "[TN]{0,1}..[ASYW][FL][CGN][LI][LI].{2,3}[RH]{2,3}[^CG]",
        "[^CG]((R[RH])|(RH))[RHQ][RH][^CG]",
        options
    )
    @test contained_complex_review_regression.matched
    @test contained_complex_review_regression.query_relationship == "Complex Overlap"
    @test contained_complex_review_regression.search_relationship == "Complex Overlap"
    @test contained_complex_review_regression.matched_pattern ==
          "[arndqehilkmfpstwyv][rh][RH][rqh][arndqehilkmfpstwyv]"
    @test contained_complex_review_regression.core_ic ≈
          0.6601809360589201 atol = 1e-12

    parent_no_longer_beats_shorter_overlap = compare(
        "[DE].{0,2}[DE].{0,2}[ADE].{0,1}[CF]\$",
        "[DE]{0,2}[DE].{0,2}[DE].{2,4}[CF].[CF]",
        options
    )
    @test parent_no_longer_beats_shorter_overlap.matched
    @test parent_no_longer_beats_shorter_overlap.query_relationship == "Complex Overlap"
    @test parent_no_longer_beats_shorter_overlap.search_relationship == "Complex Overlap"

    parent_complex_tiebreak_ignores_span_gap = compare("A.{0,1}A.{2,3}\$", "A.KA", options)
    @test parent_complex_tiebreak_ignores_span_gap.matched
    @test parent_complex_tiebreak_ignores_span_gap.query_relationship == "Complex Parent"
    @test parent_complex_tiebreak_ignores_span_gap.search_relationship ==
          "Complex Subsequence"
    @test parent_complex_tiebreak_ignores_span_gap.matched_pattern == "Aaka"

    generic_complex_shift_tiebreak = compare("[KRQ]L", "K.[KR].", options)
    @test generic_complex_shift_tiebreak.matched
    @test generic_complex_shift_tiebreak.query_relationship == "Complex Subsequence"
    @test generic_complex_shift_tiebreak.search_relationship == "Complex Parent"
    @test generic_complex_shift_tiebreak.matched_pattern == "[rqk]l"
    @test generic_complex_shift_tiebreak.core_ic ≈ 0.31663710432895764 atol = 1e-12

    reverse_generic_complex_shift_tiebreak = compare("K.[KR].", "[KRQ]L", options)
    @test reverse_generic_complex_shift_tiebreak.matched
    @test reverse_generic_complex_shift_tiebreak.query_relationship == "Complex Parent"
    @test reverse_generic_complex_shift_tiebreak.search_relationship ==
          "Complex Subsequence"
    @test reverse_generic_complex_shift_tiebreak.matched_pattern == "[rqk]l"
    @test reverse_generic_complex_shift_tiebreak.core_ic ≈ 0.3580608434035529 atol = 1e-12

    rounded_score_tiebreak = compare("H.H.{0,1}QY", "H[LM]H(([KR][^H].)|(.[^H][KR]))", options)
    @test rounded_score_tiebreak.matched
    @test rounded_score_tiebreak.query_relationship == "Complex Subsequence"
    @test rounded_score_tiebreak.search_relationship == "Complex Parent"
    @test rounded_score_tiebreak.matched_positions == 3
    @test rounded_score_tiebreak.score ≈ 1.7025165344558375 atol = 1e-12

    remapped_complex_match_tiebreak = compare(
        "[WY]..[NMR](([GQ][^I][DH][DMN])|([P][^I][DH][DMN].{0,1}[WY]))",
        "[YW].{2,4}[MNFRI].[MNF].[MNFRD]",
        options
    )
    @test remapped_complex_match_tiebreak.matched
    @test remapped_complex_match_tiebreak.query_relationship == "Complex Match"
    @test remapped_complex_match_tiebreak.search_relationship == "Complex Match"
    @test remapped_complex_match_tiebreak.matched_positions == 4
end

@testitem "alphabet specs and options surface" begin
    defaults = ComparisonOptions()
    @test defaults.alphabet isa ProteinAlphabet
    @test defaults.residue_frequencies === nothing

    rna_options = ComparisonOptions(; alphabet = RNAAlphabet())
    @test rna_options.alphabet isa RNAAlphabet

    dna_options = ComparisonOptions(; alphabet = DNAAlphabet(), min_shared_positions = 1)
    @test dna_options.alphabet isa DNAAlphabet
    @test dna_options.min_shared_positions == 1

    mismatch_options = ComparisonOptions(; mismatches = 1)
    @test mismatch_options.mismatches == 1

    custom_dna = ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict('a' => 7.0, 'c' => 2.0, 'g' => 1.0, 't' => 1.0)
    )
    @test custom_dna.residue_frequencies['A'] ≈ 7 / 11
    @test custom_dna.residue_frequencies['C'] ≈ 2 / 11
    @test custom_dna.residue_frequencies['G'] ≈ 1 / 11
    @test custom_dna.residue_frequencies['T'] ≈ 1 / 11

    built_spec = CompariMotif._build_alphabet_spec(CompariMotif.RNAAlphabet)
    cached_spec = CompariMotif._alphabet_spec(CompariMotif.RNAAlphabet())
    @test built_spec.chars == "ACGU"
    @test built_spec.index == Dict('A' => 1, 'C' => 2, 'G' => 3, 'U' => 4)
    @test built_spec.mask == UInt32(0b1111)
    @test built_spec.log_base ≈ log(4)
    @test cached_spec.chars == built_spec.chars
    @test cached_spec.index == built_spec.index
    @test cached_spec.mask == built_spec.mask
    @test cached_spec.log_base == built_spec.log_base
    @test CompariMotif._residue_frequency_vector(custom_dna, CompariMotif._alphabet_spec(DNAAlphabet())) ≈
          [7, 2, 1, 1] / 11

    huge_counts = ComparisonOptions(;
        residue_frequencies = Dict(aa => 1e308 for aa in "ACDEFGHIKLMNPQRSTVWY"),
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0
    )
    @test huge_counts.residue_frequencies['A'] ≈ 1 / 20
    huge_result = compare("A", "A", huge_counts)
    @test huge_result.matched
    @test huge_result.match_ic ≈ 1.0 atol = 1e-12
    @test huge_result.normalized_ic ≈ 1.0 atol = 1e-12
    @test isfinite(huge_result.query_information)

    shown = sprint(show, MIME"text/plain"(), ComparisonOptions(; alphabet = DNAAlphabet(), mismatches = 1))
    @test occursin("alphabet", shown)
    @test occursin("DNAAlphabet()", shown)
    @test occursin("mismatches", shown)
    @test !occursin("alphabet_index", shown)
    @test !occursin("alphabet_mask", shown)
    @test !occursin("log_base", shown)

    @test_throws ArgumentError ComparisonOptions(; min_shared_positions = 0)
    @test_throws ArgumentError ComparisonOptions(; normalized_ic_cutoff = -0.1)
    @test_throws ArgumentError ComparisonOptions(; mismatches = -1)
    @test_throws ArgumentError ComparisonOptions(; max_variants = 0)
end

@testitem "partial overlap scoring uses union information" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    spec = CompariMotif._alphabet_spec(options.alphabet)
    ugly_union = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("A[KRQ]", options), options, spec))
    ugly_query = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("A[KR]", options), options, spec))
    ugly_search = only(CompariMotif._expand_variants(
        CompariMotif._parse_motif("A[RQ]", options), options, spec))

    ugly = compare("A[KR]", "A[RQ]", options)
    @test ugly.matched
    @test ugly.query_relationship == "Complex Match"
    @test ugly.search_relationship == "Complex Match"
    @test ugly.match_ic ≈ ugly_union.information atol = 1e-12
    @test ugly.normalized_ic ≈
          ugly_union.information / min(ugly_query.information, ugly_search.information) atol = 1e-12
end

@testitem "residue frequency validation" begin
    @test_throws ArgumentError ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict('A' => 1.0, 'C' => 1.0, 'G' => 1.0)
    )
    @test_throws ArgumentError ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict(
            'A' => 1.0, 'C' => 1.0, 'G' => 1.0, 'T' => 1.0, 'U' => 1.0
        )
    )
    @test_throws ArgumentError ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict(
            'a' => 1.0, 'A' => 2.0, 'C' => 1.0, 'G' => 1.0, 'T' => 1.0
        )
    )
    @test_throws ArgumentError ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict('A' => 1.0, 'C' => NaN, 'G' => 1.0, 'T' => 1.0)
    )
    @test_throws ArgumentError ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict('A' => 1.0, 'C' => Inf, 'G' => 1.0, 'T' => 1.0)
    )
    @test_throws ArgumentError ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict('A' => 1.0, 'C' => 0.0, 'G' => 1.0, 'T' => 1.0)
    )
    @test_throws ArgumentError ComparisonOptions(;
        alphabet = DNAAlphabet(),
        residue_frequencies = Dict('A' => 1.0, 'C' => -1.0, 'G' => 1.0, 'T' => 1.0)
    )
end

@testitem "single pair relationship categories" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    exact = compare("RKLI", "RKLI", options)
    @test exact.matched
    @test exact.query_relationship == "Exact Match"
    @test exact.matched_positions == 4
    @test exact.score ≈ 4.0 atol = 1e-8

    variant = compare("RKLI", "R[KR]L[IV]", options)
    @test variant.matched
    @test variant.query_relationship == "Variant Match"
    @test variant.search_relationship == "Degenerate Match"

    degenerate = compare("R[KR]L[IV]", "RKLI", options)
    @test degenerate.matched
    @test degenerate.query_relationship == "Degenerate Match"

    complex = compare("[KR]L[IMV]", "[RKQ]L[IV]", options)
    @test complex.matched
    @test complex.query_relationship == "Complex Match"
end

@testitem "matchfix and mismatch options" begin
    base_options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    base = compare("A.", "AK", base_options)
    @test base.matched
    @test base.query_relationship == "Degenerate Match"

    qfixed_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        matchfix = :query_fixed)
    qfixed = compare("AK", "A.", qfixed_options)
    @test !qfixed.matched
    @test !compare("A[KR]L", "RL", qfixed_options).matched
    @test !compare("A.K", ".K", qfixed_options).matched
    @test !compare("A.K", "A.", qfixed_options).matched
    @test !compare("AK.A", "AK.", qfixed_options).matched
    @test compare("A[KR]L", "[KR]L", qfixed_options).matched
    @test compare(".AA", ".A", qfixed_options).matched
    @test compare("AA.", "A.", qfixed_options).matched
    @test compare("AK.A", "K.A", qfixed_options).matched
    @test compare("A.KA", "A.K", qfixed_options).matched
    @test compare(".K", "A.K", qfixed_options).matched
    qfixed_boundary = compare(".[KR]", "K.[KR].", qfixed_options)
    @test qfixed_boundary.matched
    @test qfixed_boundary.query_relationship == "Degenerate Overlap"
    @test qfixed_boundary.search_relationship == "Variant Overlap"
    @test qfixed_boundary.core_ic ≈ 0.769 atol = 1e-3

    sfixed_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        matchfix = :search_fixed)
    sfixed = compare("A.", "AK", sfixed_options)
    @test !sfixed.matched
    @test !compare("IPI", "[IPV]IAL", sfixed_options).matched
    @test !compare("RL", "A[KR]L", sfixed_options).matched
    @test !compare(".K", "A.K", sfixed_options).matched
    @test !compare("A.", "A.K", sfixed_options).matched
    @test !compare("AK.", "AK.A", sfixed_options).matched
    @test compare("[KR]L", "A[KR]L", sfixed_options).matched
    @test compare(".A", ".AA", sfixed_options).matched
    @test compare("A.", "AA.", sfixed_options).matched
    @test compare("K.A", "AK.A", sfixed_options).matched
    @test compare("A.K", "A.KA", sfixed_options).matched
    @test compare("A.K", ".K", sfixed_options).matched

    both_fixed_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        matchfix = :both_fixed)
    @test !compare("[IPV]IAL", "IPI", both_fixed_options).matched
    @test !compare("A[KR]L", "RL", both_fixed_options).matched
    @test !compare("A.K", ".K", both_fixed_options).matched
    @test !compare(".K", "A.K", both_fixed_options).matched
    @test !compare("A.K", "A.", both_fixed_options).matched
    @test !compare("A.", "A.K", both_fixed_options).matched
    @test !compare("AK.A", "AK.", both_fixed_options).matched
    @test !compare("AK.", "AK.A", both_fixed_options).matched
    @test compare("AKLI", "KLI", both_fixed_options).matched
    @test compare("A[KR]L", "[KR]L", both_fixed_options).matched
    @test compare("[IPV]IAL", "IAL", both_fixed_options).matched
    @test compare(".AA", ".A", both_fixed_options).matched
    @test compare("AA.", "A.", both_fixed_options).matched
    @test compare("AK.A", "K.A", both_fixed_options).matched
    @test compare("A.KA", "A.K", both_fixed_options).matched
    @test !compare("AKLI", "KLIA", both_fixed_options).matched

    mm0 = compare("AK", "AQ", base_options)
    @test !mm0.matched

    mm1_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        mismatches = 1)
    mm1 = compare("AK", "AQ", mm1_options)
    @test mm1.matched
    @test mm1.matched_positions == 1

    nterm_anchor_class = compare("^[KR]L", "[DE][KR]L", mm1_options)
    @test nterm_anchor_class.matched
    @test nterm_anchor_class.matched_pattern == raw"[^de][RK]L"

    anchor_shift_preferred = compare(raw"^K", raw"^AK", mm1_options)
    @test anchor_shift_preferred.matched
    @test anchor_shift_preferred.query_relationship == "Exact Subsequence"
    @test anchor_shift_preferred.search_relationship == "Exact Parent"
    @test anchor_shift_preferred.matched_pattern == raw"^[ak]"
    @test anchor_shift_preferred.core_ic ≈ 0.5 atol = 1e-3

    cterm_anchor_class = compare(raw"A[DE]$", raw"[KR]A$", mm1_options)
    @test cterm_anchor_class.matched
    @test cterm_anchor_class.matched_pattern == "A[\$de]"

    fullspan_anchor_preferred = compare(raw"A.$", raw"[KR]A$", mm1_options)
    @test fullspan_anchor_preferred.matched
    @test fullspan_anchor_preferred.query_relationship == "Degenerate Match"
    @test fullspan_anchor_preferred.search_relationship == "Variant Match"
    @test fullspan_anchor_preferred.matched_pattern == raw"[ark]a$"
    @test fullspan_anchor_preferred.core_ic ≈ 1 / 3 atol = 1e-3

    qfixed_mm1_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        mismatches = 1,
        matchfix = :query_fixed)
    @test !compare("^[KR]L", "[DE][KR]L", qfixed_mm1_options).matched
    @test compare("[DE][KR]L", "^[KR]L", qfixed_mm1_options).matched

    sfixed_mm1_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        mismatches = 1,
        matchfix = :search_fixed)
    @test compare("^[KR]L", "[DE][KR]L", sfixed_mm1_options).matched
    @test !compare("[DE][KR]L", "^[KR]L", sfixed_mm1_options).matched

    bothfixed_mm1_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        mismatches = 1,
        matchfix = :both_fixed)
    @test !compare("^[KR]L", "[DE][KR]L", bothfixed_mm1_options).matched
    @test !compare("[DE][KR]L", "^[KR]L", bothfixed_mm1_options).matched
    @test !compare(raw"^A.$", raw"A$", qfixed_mm1_options).matched
    @test !compare(raw"A$", raw"^A.$", sfixed_mm1_options).matched

    anchor_fixed_tie = compare(raw"^AK$", raw"^K", mm1_options)
    @test anchor_fixed_tie.matched
    @test anchor_fixed_tie.query_relationship == "Exact Parent"
    @test anchor_fixed_tie.search_relationship == "Exact Subsequence"
    @test anchor_fixed_tie.matched_pattern == raw"[^a]K"
end

@testitem "default option semantics" begin
    defaults = ComparisonOptions()

    @test !compare("A.", "AT", defaults).matched
    @test compare("A.", "AT", ComparisonOptions(; min_shared_positions = 1)).matched

    @test !compare("QRSTAA", "AAMNPQ", defaults).matched
    @test compare("QRSTAA", "AAMNPQ", ComparisonOptions(; normalized_ic_cutoff = 0.0)).matched

    @test compare("AKL", "A.L", defaults).matched
    @test !compare("AKL", "A.L", ComparisonOptions(; matchfix = :query_fixed)).matched

    @test !compare("AKL", "AQL", defaults).matched
    @test compare("AKL", "AQL", ComparisonOptions(; mismatches = 1)).matched

    overlap_default = compare("A[KR]", "A[RQ]", defaults)
    @test overlap_default.matched
    @test overlap_default.query_relationship == "Complex Match"
    boundary_overlap = compare(
        ".[KR]", "K.[KR].", ComparisonOptions(;
            min_shared_positions = 1,
            normalized_ic_cutoff = 0.0
        ))
    @test boundary_overlap.matched
    @test boundary_overlap.query_relationship == "Degenerate Overlap"
    @test boundary_overlap.search_relationship == "Variant Overlap"
    @test boundary_overlap.core_ic ≈ 0.769 atol = 1e-3
    reverse_boundary = compare(
        "K.[KR].", ".[KR]", ComparisonOptions(;
            min_shared_positions = 1,
            normalized_ic_cutoff = 0.0
        ))
    @test reverse_boundary.matched
    @test reverse_boundary.query_relationship == "Exact Parent"
    @test reverse_boundary.search_relationship == "Exact Subsequence"
    @test !compare(
        "A[KR]",
        "A[RQ]",
        ComparisonOptions(; allow_ambiguous_overlap = false)
    ).matched
end

@testitem "matchfix validation" begin
    @test ComparisonOptions(; matchfix = :none).matchfix == :none
    @test ComparisonOptions(; matchfix = :query_fixed).matchfix == :query_fixed
    @test ComparisonOptions(; matchfix = :search_fixed).matchfix == :search_fixed
    @test ComparisonOptions(; matchfix = :both_fixed).matchfix == :both_fixed
    @test ComparisonOptions(
        ProteinAlphabet(),
        nothing,
        2,
        0.5,
        :both_fixed,
        0,
        true,
        10_000
    ).matchfix == :both_fixed
    @test ComparisonOptions(
        ProteinAlphabet(),
        nothing,
        2,
        0,
        :none,
        0,
        true,
        10_000
    ).normalized_ic_cutoff == 0.0

    integer_frequencies = Dict(aa => 1 for aa in collect("ACDEFGHIKLMNPQRSTVWY"))
    positional_frequencies = ComparisonOptions(
        ProteinAlphabet(),
        integer_frequencies,
        2,
        0.5,
        :none,
        0,
        true,
        10_000
    ).residue_frequencies
    @test positional_frequencies isa Dict{Char, Float64}
    @test positional_frequencies !== nothing
    @test isapprox(sum(values(positional_frequencies)), 1.0)
    @test all(values(positional_frequencies) .≈ 0.05)

    expected_message = "`matchfix` must be one of :none, :query_fixed, :search_fixed, :both_fixed."

    invalid_symbol = try
        ComparisonOptions(; matchfix = :query)
        nothing
    catch err
        err
    end
    @test invalid_symbol isa ArgumentError
    @test occursin(expected_message, sprint(showerror, invalid_symbol))

    invalid_string = try
        ComparisonOptions(; matchfix = "query_fixed")
        nothing
    catch err
        err
    end
    @test invalid_string isa TypeError

    invalid_positional = try
        ComparisonOptions(ProteinAlphabet(), nothing, 2, 0.5, :query, 0, true, 10_000)
        nothing
    catch err
        err
    end
    @test invalid_positional isa ArgumentError
    @test occursin(expected_message, sprint(showerror, invalid_positional))

    invalid_runtime = try
        CompariMotif._query_fixed_required(:query)
        nothing
    catch err
        err
    end
    @test invalid_runtime isa ArgumentError
    @test occursin(expected_message, sprint(showerror, invalid_runtime))

    invalid_search_runtime = try
        CompariMotif._search_fixed_required(:query)
        nothing
    catch err
        err
    end
    @test invalid_search_runtime isa ArgumentError
    @test occursin(expected_message, sprint(showerror, invalid_search_runtime))
end

@testitem "candidate tie-break ordering" begin
    candidate = (; matched_positions,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 1,
        dual_wildcard_positions = 0,
        exact_positions = 0,
        leading_exact_positions = 0,
        rel_type = CompariMotif._REL_EXACT,
        rel_length = CompariMotif._LEN_MATCH,
        query_rel_type = rel_type,
        query_rel_length = rel_length,
        search_rel_type = rel_type,
        search_rel_length = rel_length,
        query_normalized = "Q",
        search_normalized = "S",
        query_information = 1.0,
        search_information = 1.0) -> CompariMotif._Candidate(
        CompariMotif._MotifVariant(CompariMotif._Position[], query_normalized, query_information),
        CompariMotif._MotifVariant(CompariMotif._Position[], search_normalized, search_information),
        query_rel_type,
        query_rel_length,
        search_rel_type,
        search_rel_length,
        "",
        matched_positions,
        match_ic,
        1.0,
        core_ic,
        score,
        overlap_length,
        dual_wildcard_positions,
        exact_positions,
        leading_exact_positions
    )

    higher_match_ic = candidate(; matched_positions = 2, match_ic = 2.0, score = 1.0)
    lower_match_ic = candidate(; matched_positions = 3, match_ic = 1.0, score = 3.0)
    @test CompariMotif._is_better(higher_match_ic, lower_match_ic)

    more_positions = candidate(; matched_positions = 3, match_ic = 1.0, score = 1.0)
    fewer_positions = candidate(; matched_positions = 2, match_ic = 1.0, score = 2.0)
    @test CompariMotif._is_better(more_positions, fewer_positions)

    higher_score = candidate(; matched_positions = 3, match_ic = 1.0, score = 2.0)
    lower_score = candidate(; matched_positions = 3, match_ic = 1.0, score = 1.0)
    @test CompariMotif._is_better(higher_score, lower_score)

    tied_a = candidate(; matched_positions = 3, match_ic = 1.0, score = 1.0)
    tied_b = candidate(; matched_positions = 3, match_ic = 1.0, score = 1.0)
    @test !CompariMotif._is_better(tied_a, tied_b)

    fewer_dual_wildcards = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        dual_wildcard_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_MATCH)
    more_dual_wildcards = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        dual_wildcard_positions = 1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_MATCH)
    @test CompariMotif._is_better(fewer_dual_wildcards, more_dual_wildcards)

    higher_core = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT)
    lower_core = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP)
    @test !CompariMotif._is_better(higher_core, lower_core)
    @test !CompariMotif._is_better(lower_core, higher_core)

    clean_overlap = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 3,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP,
        query_normalized = "ABCDE")
    medium_parent = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 4,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT,
        query_normalized = "ABCDE")
    cleaner_subsequence = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 5,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_SUBSEQUENCE,
        query_normalized = "ABCDE")
    longer_overlap = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 6,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP,
        query_normalized = "ABCDE")
    padding_heavier_subsequence = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 5,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_SUBSEQUENCE,
        query_normalized = "ABCDE")
    @test !CompariMotif._is_better(cleaner_subsequence, longer_overlap)
    @test !CompariMotif._is_better(longer_overlap, cleaner_subsequence)
    @test !CompariMotif._is_better(clean_overlap, padding_heavier_subsequence)
    @test !CompariMotif._is_better(padding_heavier_subsequence, clean_overlap)
    @test !CompariMotif._is_better(clean_overlap, medium_parent)
    @test !CompariMotif._is_better(medium_parent, clean_overlap)

    exacter_subsequence = candidate(;
        matched_positions = 4,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 6,
        exact_positions = 3,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_SUBSEQUENCE,
        query_normalized = "ABCDE")
    looser_overlap = candidate(;
        matched_positions = 4,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 5,
        exact_positions = 1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP,
        query_normalized = "ABCDE")
    @test !CompariMotif._is_better(exacter_subsequence, looser_overlap)
    @test !CompariMotif._is_better(looser_overlap, exacter_subsequence)

    parent_with_longer_exact_prefix = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 4,
        exact_positions = 1,
        leading_exact_positions = 1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT)
    overlap_without_exact_prefix = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 3,
        exact_positions = 1,
        leading_exact_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP)
    @test CompariMotif._is_better(parent_with_longer_exact_prefix, overlap_without_exact_prefix)
    @test !CompariMotif._is_better(
        overlap_without_exact_prefix,
        parent_with_longer_exact_prefix
    )

    parent_without_prefix_advantage = candidate(;
        matched_positions = 2,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 7,
        exact_positions = 2,
        leading_exact_positions = 1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT)
    overlap_with_same_prefix = candidate(;
        matched_positions = 2,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.1,
        overlap_length = 6,
        exact_positions = 2,
        leading_exact_positions = 1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP)
    @test !CompariMotif._is_better(overlap_with_same_prefix, parent_without_prefix_advantage)
    @test !CompariMotif._is_better(
        parent_without_prefix_advantage,
        overlap_with_same_prefix
    )

    single_position_subsequence = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 5,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_SUBSEQUENCE)
    single_position_overlap = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.03,
        overlap_length = 4,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP)
    @test !CompariMotif._is_better(single_position_overlap, single_position_subsequence)
    @test !CompariMotif._is_better(single_position_subsequence, single_position_overlap)

    same_span_lower_dual_subsequence = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 3,
        dual_wildcard_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_SUBSEQUENCE)
    same_span_higher_dual_overlap = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.2,
        overlap_length = 3,
        dual_wildcard_positions = 1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP)
    @test !CompariMotif._is_better(
        same_span_lower_dual_subsequence,
        same_span_higher_dual_overlap
    )
    @test !CompariMotif._is_better(
        same_span_higher_dual_overlap,
        same_span_lower_dual_subsequence
    )

    same_variant_higher_core = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.2,
        overlap_length = 3,
        dual_wildcard_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT,
        search_rel_length = CompariMotif._LEN_MATCH,
        query_normalized = "Qsame",
        search_normalized = "Ssame")
    same_variant_lower_core = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 3,
        dual_wildcard_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT,
        search_rel_length = CompariMotif._LEN_PARENT,
        query_normalized = "Qsame",
        search_normalized = "Ssame")
    @test CompariMotif._is_better(same_variant_higher_core, same_variant_lower_core)
    @test !CompariMotif._is_better(same_variant_lower_core, same_variant_higher_core)

    generic_complex_earlier = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 0.317,
        overlap_length = 2,
        dual_wildcard_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_SUBSEQUENCE)
    generic_complex_later_higher_core = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 0.358,
        overlap_length = 2,
        dual_wildcard_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_SUBSEQUENCE)
    @test !CompariMotif._is_better(generic_complex_later_higher_core, generic_complex_earlier)
    @test !CompariMotif._is_better(generic_complex_earlier, generic_complex_later_higher_core)

    shorter_span_different_signature = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 3,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT,
        search_rel_length = CompariMotif._LEN_OVERLAP)
    longer_span_different_signature = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 4,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT,
        search_rel_length = CompariMotif._LEN_PARENT)
    @test CompariMotif._is_better(
        shorter_span_different_signature,
        longer_span_different_signature
    )
    @test !CompariMotif._is_better(
        longer_span_different_signature,
        shorter_span_different_signature
    )

    cross_variant_same_signature_earlier = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        dual_wildcard_positions = 1,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP,
        query_normalized = "Q1")
    cross_variant_same_signature_later = candidate(;
        matched_positions = 1,
        match_ic = 1.0,
        score = 1.0,
        dual_wildcard_positions = 0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_OVERLAP,
        query_normalized = "Q2")
    @test !CompariMotif._is_better(
        cross_variant_same_signature_later,
        cross_variant_same_signature_earlier
    )
    @test !CompariMotif._is_better(
        cross_variant_same_signature_earlier,
        cross_variant_same_signature_later
    )

    shorter_contained_overlap = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 4,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT)
    longer_contained_overlap = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        core_ic = 1.0,
        overlap_length = 5,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT)
    @test !CompariMotif._is_better(shorter_contained_overlap, longer_contained_overlap)

    rounded_higher = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0 + 5e-13,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_MATCH)
    rounded_lower = candidate(;
        matched_positions = 3,
        match_ic = 1.0,
        score = 1.0,
        rel_type = CompariMotif._REL_COMPLEX,
        rel_length = CompariMotif._LEN_PARENT)
    @test !CompariMotif._is_better(rounded_higher, rounded_lower)

    exact_full_lower_information = candidate(;
        matched_positions = 2,
        match_ic = 2.0,
        score = 2.0,
        rel_type = CompariMotif._REL_EXACT,
        rel_length = CompariMotif._LEN_MATCH,
        query_information = 2.0,
        search_information = 2.0)
    exact_parent_higher_information = candidate(;
        matched_positions = 2,
        match_ic = 2.0,
        score = 2.0,
        rel_type = CompariMotif._REL_EXACT,
        rel_length = CompariMotif._LEN_PARENT,
        query_information = 2.77,
        search_information = 2.0)
    @test CompariMotif._is_better(
        exact_full_lower_information,
        exact_parent_higher_information
    )
    @test !CompariMotif._is_better(
        exact_parent_higher_information,
        exact_full_lower_information
    )

    exact_full_equal_information = candidate(;
        matched_positions = 2,
        match_ic = 2.0,
        score = 2.0,
        rel_type = CompariMotif._REL_EXACT,
        rel_length = CompariMotif._LEN_MATCH,
        query_information = 2.0,
        search_information = 2.0)
    exact_parent_equal_information = candidate(;
        matched_positions = 2,
        match_ic = 2.0,
        score = 2.0,
        rel_type = CompariMotif._REL_EXACT,
        rel_length = CompariMotif._LEN_PARENT,
        query_information = 2.0,
        search_information = 2.0)
    @test !CompariMotif._is_better(
        exact_full_equal_information,
        exact_parent_equal_information
    )
end

@testitem "normalized-ic cutoff keeps boundary hits" begin
    options = ComparisonOptions()
    residue = CompariMotif._Position(CompariMotif._RESIDUE, CompariMotif.ResidueMask(1))
    query_variant = CompariMotif._MotifVariant(fill(residue, 4), "Q", 3.537243573680482)
    search_variant = CompariMotif._MotifVariant(fill(residue, 8), "S", 5.7686217868402405)
    acc = CompariMotif._AlignmentAccumulator(;
        matched_positions = 2,
        match_ic = 1.7686217868402407,
        core_ic_denominator = 2.0
    )
    candidate = CompariMotif._build_alignment_candidate(
        query_variant,
        search_variant,
        2,
        acc,
        options
    )

    @test candidate !== nothing
    @test candidate.match_ic ≈ 1.7686217868402407 atol = 1e-12
    @test candidate.normalized_ic ≈ 0.5 atol = 1e-12

    denom = min(query_variant.information, search_variant.information)
    near_cutoff_acc = CompariMotif._AlignmentAccumulator(;
        matched_positions = 2,
        match_ic = denom * (options.normalized_ic_cutoff - 5e-13),
        core_ic_denominator = 2.0
    )
    near_cutoff = CompariMotif._build_alignment_candidate(
        query_variant,
        search_variant,
        2,
        near_cutoff_acc,
        options
    )
    @test near_cutoff !== nothing
    @test near_cutoff.normalized_ic ≈ (options.normalized_ic_cutoff - 5e-13) atol = 1e-15

    below_cutoff_acc = CompariMotif._AlignmentAccumulator(;
        matched_positions = 2,
        match_ic = denom * (options.normalized_ic_cutoff - 2e-12),
        core_ic_denominator = 2.0
    )
    below_cutoff = CompariMotif._build_alignment_candidate(
        query_variant,
        search_variant,
        2,
        below_cutoff_acc,
        options
    )
    @test isnothing(below_cutoff)
end

@testitem "matrix APIs" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    motifs = ["RKLI", "R[KR]L[IV]", "[KR].L.[FYLIMVP]"]
    matrix_self = compare(motifs, options)
    @test size(matrix_self) == (3, 3)

    matrix_db = compare(motifs, ["R.LE", "RKL"], options)
    @test size(matrix_db) == (3, 2)
    @test matrix_db[3, 1].matched
    @test matrix_db[3, 1].query_relationship == "Degenerate Parent"

    single_vs_many = compare("RKLI", ["RKLI", "R[KR]L[IV]"], options)
    wrapped_single_vs_many = compare(["RKLI"], ["RKLI", "R[KR]L[IV]"], options)
    @test size(single_vs_many) == (1, 2)
    @test to_column_table(single_vs_many) == to_column_table(wrapped_single_vs_many)
end

@testitem "options overloads and abstract string interfaces" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    pair = compare("RKLI", "R[KR]L[IV]", options)
    @test pair.matched
    @test pair.query_relationship == "Variant Match"

    split_motifs = split("RKLI,R[KR]L[IV],[KR].L.[FYLIMVP]", ',')
    @test split_motifs isa Vector{<:AbstractString}
    matrix = compare(split_motifs, options)
    @test size(matrix) == (3, 3)
    @test matrix[2, 1].query_relationship == "Degenerate Match"

    @test_throws MethodError compare(split_motifs, "RKLI", options)
end

@testitem "max_variants overflow guard" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    repeated_choices = (8 * sizeof(Int)) - 1
    overflow_motif = repeat("x{1,2}", repeated_choices)
    parsed = CompariMotif._parse_motif(overflow_motif, options)
    nvariants = CompariMotif._variant_count(parsed.tokens)

    @test nvariants == big(2)^repeated_choices
    @test nvariants > big(typemax(Int))

    large_but_manageable = repeat("x{1,10}", 6)  # 10^6 variants
    strict_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        max_variants = 100_000
    )
    @test_throws ArgumentError compare(large_but_manageable, "A", strict_options)

    strict_alternation_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        max_variants = 100
    )
    @test_throws ArgumentError CompariMotif._parse_motif(
        repeat("(A|C)", 7),
        strict_alternation_options
    )
end

@testitem "to_column_table overloads" begin
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    pair = compare("RKLI", "R[KR]L[IV]", options)

    single = to_column_table(pair)
    @test keys(single) == (
        :query, :search, :normalized_query, :normalized_search, :matched,
        :query_relationship, :search_relationship, :matched_pattern,
        :matched_positions, :match_ic, :normalized_ic, :core_ic, :score,
        :query_information, :search_information
    )
    @test length(single.query) == 1
    @test single.query_relationship[1] == "Variant Match"

    vec_table = to_column_table([pair, pair])
    @test keys(vec_table)[1] == :result_index
    @test vec_table.result_index == [1, 2]
    @test vec_table.search_relationship == ["Degenerate Match", "Degenerate Match"]

    matrix = compare(["RKLI", "R[KR]L[IV]"], options)
    mat_table = to_column_table(matrix)
    @test keys(mat_table)[1:2] == (:query_index, :search_index)
    @test length(mat_table.query_index) == 4
    @test mat_table.query_index == [1, 1, 2, 2]
    @test mat_table.search_index == [1, 2, 1, 2]

    single_vs_many = to_column_table(compare("RKLI", ["RKLI", "R[KR]L[IV]"], options))
    @test single_vs_many.query_index == [1, 1]
    @test single_vs_many.search_index == [1, 2]
end

@testitem "to_column_table DataFrame conversion" begin
    import DataFrames

    motifs = ["RKLI", "R[KR]L[IV]"]
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    df = DataFrames.DataFrame(to_column_table(compare(motifs, options)))

    @test size(df) == (4, 17)
    @test names(df)[1:5] ==
          ["query_index", "search_index", "query", "search", "normalized_query"]
    @test df.query_relationship == [
        "Exact Match",
        "Variant Match",
        "Degenerate Match",
        "Exact Match"
    ]
end

@testitem "to_column_table CSV export" begin
    import CSV

    motifs = ["RKLI", "R[KR]L[IV]"]
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    table = to_column_table(compare(motifs, options))

    mktemp() do path, io
        close(io)
        CSV.write(path, table)
        @test isfile(path)
        @test filesize(path) > 0

        rows = collect(CSV.File(path))
        @test length(rows) == 4
        @test rows[2].query_relationship == "Variant Match"
    end
end

@testitem "RNA alphabet mode" begin
    rna_options = ComparisonOptions(;
        alphabet = RNAAlphabet(), min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test CompariMotif._normalize_motif("A[CU]x"; alphabet = RNAAlphabet()) == "A[CU]."
    @test compare("AUG", "AUG", rna_options).matched
    @test_throws ArgumentError CompariMotif._normalize_motif("ATG"; alphabet = RNAAlphabet())
end
