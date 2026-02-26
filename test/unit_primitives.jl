using TestItems

@testitem "normalize motif syntax" begin
    using Test
    using CompariMotif

    @test normalize_motif("r[kR].{0,1}l") == "R[KR]x{0,1}L"
    @test normalize_motif("A[^P]x") == "A[ACDEFGHIKLMNQRSTVWY]x"
    @test normalize_motif("x(1,2)") == "x{1,2}"
end

@testitem "single pair relationship categories" begin
    using Test
    using CompariMotif

    exact = compare("RKLI", "RKLI")
    @test exact.matched
    @test exact.query_relationship == "Exact Match"
    @test exact.query_relationship_code == "e-m"
    @test exact.matched_positions == 4
    @test exact.score ≈ 4.0 atol = 1e-8

    variant = compare("RKLI", "R[KR]L[IV]"; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test variant.matched
    @test variant.query_relationship == "Variant Match"
    @test variant.search_relationship == "Degenerate Match"
    @test variant.query_relationship_code == "v-m"

    degenerate = compare("R[KR]L[IV]", "RKLI"; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test degenerate.matched
    @test degenerate.query_relationship == "Degenerate Match"
    @test degenerate.query_relationship_code == "d-m"

    complex = compare("[KR]L[IMV]", "[RKQ]L[IV]"; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test complex.matched
    @test complex.query_relationship == "Complex Match"
    @test complex.query_relationship_code == "c-m"
end

@testitem "matchfix and mismatch options" begin
    using Test
    using CompariMotif

    base = compare("A.", "AK"; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test base.matched
    @test base.query_relationship == "Degenerate Match"

    qfixed = compare(
        "AK", "A."; min_shared_positions = 1, normalized_ic_cutoff = 0.0, matchfix = 1)
    @test !qfixed.matched

    sfixed = compare(
        "A.", "AK"; min_shared_positions = 1, normalized_ic_cutoff = 0.0, matchfix = 2)
    @test !sfixed.matched

    mm0 = compare("AK", "AQ"; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test !mm0.matched

    mm1 = compare(
        "AK", "AQ"; min_shared_positions = 1, normalized_ic_cutoff = 0.0, mismatches = 1)
    @test mm1.matched
    @test mm1.matched_positions == 1
end

@testitem "matrix APIs" begin
    using Test
    using CompariMotif

    motifs = ["RKLI", "R[KR]L[IV]", "[KR]xLx[FYLIMVP]"]
    matrix_self = compare(motifs; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test size(matrix_self) == (3, 3)

    matrix_db = compare(motifs, ["RxLE", "RKL"]; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    @test size(matrix_db) == (3, 2)
    @test matrix_db[3, 1].matched
    @test matrix_db[3, 1].query_relationship == "Degenerate Parent"
end

@testitem "max_variants overflow guard" begin
    using Test
    using CompariMotif

    options = CompariMotif._build_options(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    repeated_choices = (8 * sizeof(Int)) - 1
    overflow_motif = repeat("x{1,2}", repeated_choices)
    parsed = CompariMotif._parse_motif(overflow_motif, options)
    nvariants = CompariMotif._variant_count(parsed.tokens)

    @test nvariants == big(2)^repeated_choices
    @test nvariants > big(typemax(Int))

    large_but_manageable = repeat("x{1,10}", 6)  # 10^6 variants
    @test_throws ArgumentError compare(
        large_but_manageable,
        "A";
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        max_variants = 100_000
    )
end
