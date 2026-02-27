using TestItems

@testitem "normalize motif syntax" begin
    using Test
    using CompariMotif

    @test normalize_motif("r[kR].{0,1}l") == "R[RK]x{0,1}L"
    @test normalize_motif("A[^P]x") == "A[ARNDCQEGHILKMFSTWYV]x"
    @test normalize_motif("x(1,2)") == "x{1,2}"
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

@testitem "RNA alphabet mode" begin
    using Test
    using CompariMotif

    rna_options = ComparisonOptions(; alphabet = :rna, min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test normalize_motif("A[CU]x"; alphabet = :rna) == "A[CU]x"
    @test compare("AUG", "AUG", rna_options).matched
    @test_throws ArgumentError normalize_motif("ATG"; alphabet = :rna)
end
