using TestItems

@testitem "normalize motif syntax" begin
    using Test
    using CompariMotif

    @test normalize_motif("r[kR].{0,1}l") == "R[RK]x{0,1}L"
    @test normalize_motif("A[^P]x") == "A[ARNDCQEGHILKMFSTWYV]x"
    @test_throws ArgumentError normalize_motif("x(1,2)")
end

@testitem "alternation syntax" begin
    using Test
    using CompariMotif

    @test normalize_motif("(rkli)") == "RKLI"
    @test normalize_motif("(rkli)|(r[kr]l[iv])") == "(RKLI)|(R[RK]L[IV])"
    @test normalize_motif("RKLI|R[KR]L[IV]") == "(RKLI)|(R[RK]L[IV])"
    @test normalize_motif("A(K|Q)LI") == "(AKLI)|(AQLI)"
    @test normalize_motif("(K|Q)") == "(K)|(Q)"
    @test normalize_motif("R(KL)I") == "RKLI"

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    result = compare("(RKLI)|(AQLI)", "AQLI", options)
    @test result.matched
    @test result.query_relationship == "Exact Match"
    @test result.normalized_query == "(RKLI)|(AQLI)"

    class_equivalence = compare("(K|Q)", "[KQ]", options)
    @test class_equivalence.matched
    @test class_equivalence.query_relationship == "Variant Match"
    @test class_equivalence.search_relationship == "Degenerate Match"

    redundant_grouping = compare("R(KL)I", "RKLI", options)
    @test redundant_grouping.matched
    @test redundant_grouping.query_relationship == "Exact Match"
    @test redundant_grouping.normalized_query == "RKLI"
end

@testitem "wildcard token equivalence" begin
    using Test
    using CompariMotif

    @test normalize_motif("A.Xx"; alphabet = ProteinAlphabet()) == "Axxx"
    @test normalize_motif("A.Xx"; alphabet = DNAAlphabet()) == "Axxx"
    @test normalize_motif("A.Xx"; alphabet = RNAAlphabet()) == "Axxx"
end

@testitem "position information content" begin
    using Test
    using CompariMotif

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

@testitem "non-uniform residue frequencies" begin
    using Test
    using CompariMotif

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
    using Test
    using CompariMotif

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
end

@testitem "alphabet specs and options surface" begin
    using Test
    using CompariMotif

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
end

@testitem "partial overlap scoring uses union information" begin
    using Test
    using CompariMotif

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
    using Test
    using CompariMotif

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
    using Test
    using CompariMotif

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
    using Test
    using CompariMotif

    base_options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    base = compare("A.", "AK", base_options)
    @test base.matched
    @test base.query_relationship == "Degenerate Match"

    qfixed_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        matchfix = MatchFixQueryFixed)
    qfixed = compare("AK", "A.", qfixed_options)
    @test !qfixed.matched

    sfixed_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        matchfix = MatchFixSearchFixed)
    sfixed = compare("A.", "AK", sfixed_options)
    @test !sfixed.matched

    mm0 = compare("AK", "AQ", base_options)
    @test !mm0.matched

    mm1_options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        mismatches = 1)
    mm1 = compare("AK", "AQ", mm1_options)
    @test mm1.matched
    @test mm1.matched_positions == 1
end

@testitem "default option semantics" begin
    using Test
    using CompariMotif

    defaults = ComparisonOptions()

    @test !compare("A.", "AT", defaults).matched
    @test compare("A.", "AT", ComparisonOptions(; min_shared_positions = 1)).matched

    @test !compare("QRSTAA", "AAMNPQ", defaults).matched
    @test compare("QRSTAA", "AAMNPQ", ComparisonOptions(; normalized_ic_cutoff = 0.0)).matched

    @test compare("AKL", "A.L", defaults).matched
    @test !compare("AKL", "A.L", ComparisonOptions(; matchfix = MatchFixQueryFixed)).matched

    @test !compare("AKL", "AQL", defaults).matched
    @test compare("AKL", "AQL", ComparisonOptions(; mismatches = 1)).matched

    overlap_default = compare("A[KR]", "A[RQ]", defaults)
    @test overlap_default.matched
    @test overlap_default.query_relationship == "Complex Match"
    @test !compare(
        "A[KR]",
        "A[RQ]",
        ComparisonOptions(; allow_ambiguous_overlap = false)
    ).matched
end

@testitem "matrix APIs" begin
    using Test
    using CompariMotif

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    motifs = ["RKLI", "R[KR]L[IV]", "[KR]xLx[FYLIMVP]"]
    matrix_self = compare(motifs, options)
    @test size(matrix_self) == (3, 3)

    matrix_db = compare(motifs, ["RxLE", "RKL"], options)
    @test size(matrix_db) == (3, 2)
    @test matrix_db[3, 1].matched
    @test matrix_db[3, 1].query_relationship == "Degenerate Parent"
end

@testitem "options overloads and abstract string interfaces" begin
    using Test
    using CompariMotif

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    pair = compare("RKLI", "R[KR]L[IV]", options)
    @test pair.matched
    @test pair.query_relationship == "Variant Match"

    split_motifs = split("RKLI,R[KR]L[IV],[KR]xLx[FYLIMVP]", ',')
    @test split_motifs isa Vector{<:AbstractString}
    matrix = compare(split_motifs, options)
    @test size(matrix) == (3, 3)
    @test matrix[2, 1].query_relationship == "Degenerate Match"
end

@testitem "max_variants overflow guard" begin
    using Test
    using CompariMotif

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
end

@testitem "to_column_table overloads" begin
    using Test
    using CompariMotif

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
end

@testitem "to_column_table DataFrame conversion" begin
    using Test
    using DataFrames
    using CompariMotif

    motifs = ["RKLI", "R[KR]L[IV]"]
    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    df = DataFrame(to_column_table(compare(motifs, options)))

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
    using Test
    using CSV
    using CompariMotif

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
    using Test
    using CompariMotif

    rna_options = ComparisonOptions(;
        alphabet = RNAAlphabet(), min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test normalize_motif("A[CU]x"; alphabet = RNAAlphabet()) == "A[CU]x"
    @test compare("AUG", "AUG", rna_options).matched
    @test_throws ArgumentError normalize_motif("ATG"; alphabet = RNAAlphabet())
end
