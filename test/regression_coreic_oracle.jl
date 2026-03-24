@testitem "oracle CoreIC parity fixtures" begin
    using CompariMotif: CompariMotif, ComparisonOptions, ProteinAlphabet, compare

    function load_fixture_motifs(path::String)
        motifs = String[]
        for line in eachline(path)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, '#') && continue
            cols = split(stripped, '\t')
            length(cols) == 1 || error("Invalid motif fixture line: $line")
            push!(motifs, cols[1])
        end
        return motifs
    end

    function load_oracle_rows(path::String)
        header = String[]
        rows = Dict{String, String}[]
        for line in eachline(path)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, '#') && continue
            cols = split(stripped, '\t')
            if isempty(header)
                header = cols
                continue
            end
            length(cols) == length(header) || error("Malformed oracle fixture row: $line")
            row = Dict{String, String}()
            for (key, value) in zip(header, cols)
                row[key] = value
            end
            push!(rows, row)
        end
        return rows
    end

    function canonicalize_match_pattern(pattern::AbstractString)
        spec = CompariMotif._alphabet_spec(ProteinAlphabet())
        io = IOBuffer()
        i = firstindex(pattern)
        while i <= lastindex(pattern)
            if pattern[i] == '['
                close_idx = findnext(==(']'), pattern, nextind(pattern, i))
                close_idx === nothing && error("Malformed match pattern: $pattern")
                body = pattern[nextind(pattern, i):prevind(pattern, close_idx)]
                anchors = String(filter(c -> c == '^' || c == '$', body))
                residues = String(filter(c -> c != '^' && c != '$', body))
                normalized_residues = if isempty(residues)
                    ""
                else
                    mask = CompariMotif._class_mask(uppercase(residues), spec)
                    CompariMotif._mask_to_symbol(
                        mask,
                        spec;
                        as_lowercase = all(islowercase, residues)
                    )
                end
                print(io, '[', anchors, normalized_residues, ']')
                i = nextind(pattern, close_idx)
                continue
            end
            print(io, pattern[i])
            i = nextind(pattern, i)
        end
        return String(take!(io))
    end

    oracle_to_package(rel::String) = replace(rel, "Ugly" => "Complex")

    root = dirname(@__DIR__)
    motifs_path = joinpath(root, "data", "fixtures", "coreic_probe_motifs.tsv")
    patterns = load_fixture_motifs(motifs_path)
    motif_to_index = Dict(pattern => idx for (idx, pattern) in pairs(patterns))
    @test length(motif_to_index) == length(patterns)

    fixture_specs = [
        (
            oracle_path = joinpath(
                root,
                "data",
                "fixtures",
                "oracle_coreic_probe_normalized.tsv"
            ),
            options = ComparisonOptions(;
                min_shared_positions = 1,
                normalized_ic_cutoff = 0.0
            )
        ),
        (
            oracle_path = joinpath(
                root,
                "data",
                "fixtures",
                "oracle_coreic_mismatch_probe_normalized.tsv"
            ),
            options = ComparisonOptions(;
                min_shared_positions = 1,
                normalized_ic_cutoff = 0.0,
                mismatches = 1
            )
        )
    ]

    for spec in fixture_specs
        results = compare(patterns, patterns, spec.options)
        oracle_rows = load_oracle_rows(spec.oracle_path)

        expected_pairs = Set{Tuple{Int, Int}}()
        for row in oracle_rows
            qi = motif_to_index[row["Motif1"]]
            si = motif_to_index[row["Motif2"]]
            push!(expected_pairs, (qi, si))

            result = results[qi, si]
            @test result.matched
            @test result.query_relationship == oracle_to_package(row["Sim1"])
            @test result.search_relationship == oracle_to_package(row["Sim2"])
            @test canonicalize_match_pattern(result.matched_pattern) ==
                  canonicalize_match_pattern(row["Match"])
            @test result.matched_positions == parse(Int, row["MatchPos"])
            @test result.match_ic ≈ parse(Float64, row["MatchIC"]) atol = 1e-3
            @test result.normalized_ic ≈ parse(Float64, row["NormIC"]) atol = 1e-3
            @test result.core_ic ≈ parse(Float64, row["CoreIC"]) atol = 1e-3
            @test result.score ≈ parse(Float64, row["Score"]) atol = 1e-3
            @test result.query_information ≈ parse(Float64, row["Info1"]) atol = 1e-2
            @test result.search_information ≈ parse(Float64, row["Info2"]) atol = 1e-2
        end
        @test length(expected_pairs) == length(oracle_rows)

        observed_pairs = Set{Tuple{Int, Int}}()
        for i in axes(results, 1), j in axes(results, 2)

            if i != j && results[i, j].matched
                push!(observed_pairs, (i, j))
            end
        end
        @test observed_pairs == expected_pairs
    end
end

@testitem "exact self-match CoreIC probe" begin
    using CompariMotif: ComparisonOptions, compare

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    result = compare(raw"A..K", raw"A..K", options)

    # These values were obtained from a black-box CompariMotif V3.14.1 oracle run
    # against identical A..K motifs in separate input/search files.
    @test result.matched
    @test result.query_relationship == "Exact Match"
    @test result.search_relationship == "Exact Match"
    @test result.matched_pattern == raw"A..K"
    @test result.matched_positions == 2
    @test result.match_ic ≈ 2.000 atol = 1e-3
    @test result.normalized_ic ≈ 1.000 atol = 1e-3
    @test result.core_ic ≈ 1.000 atol = 1e-3
    @test result.score ≈ 2.000 atol = 1e-3
    @test result.query_information ≈ 2.00 atol = 1e-2
    @test result.search_information ≈ 2.00 atol = 1e-2
end

@testitem "mismatch-tolerant anchor CoreIC probes" begin
    using CompariMotif: CompariMotif, ComparisonOptions, ProteinAlphabet, compare

    function canonicalize_match_pattern(pattern::AbstractString)
        spec = CompariMotif._alphabet_spec(ProteinAlphabet())
        io = IOBuffer()
        i = firstindex(pattern)
        while i <= lastindex(pattern)
            if pattern[i] == '['
                close_idx = findnext(==(']'), pattern, nextind(pattern, i))
                close_idx === nothing && error("Malformed match pattern: $pattern")
                body = pattern[nextind(pattern, i):prevind(pattern, close_idx)]
                anchors = String(filter(c -> c == '^' || c == '$', body))
                residues = String(filter(c -> c != '^' && c != '$', body))
                normalized_residues = if isempty(residues)
                    ""
                else
                    mask = CompariMotif._class_mask(uppercase(residues), spec)
                    CompariMotif._mask_to_symbol(
                        mask,
                        spec;
                        as_lowercase = all(islowercase, residues)
                    )
                end
                print(io, '[', anchors, normalized_residues, ']')
                i = nextind(pattern, close_idx)
                continue
            end
            print(io, pattern[i])
            i = nextind(pattern, i)
        end
        return String(take!(io))
    end

    options = ComparisonOptions(;
        min_shared_positions = 1,
        normalized_ic_cutoff = 0.0,
        mismatches = 1
    )
    default_options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    oracle_cases = [
        (
            query = raw"A..K",
            search = raw"^AK$",
            query_relationship = "Degenerate Overlap",
            search_relationship = "Variant Overlap",
            matched_pattern = raw"^aK",
            matched_positions = 1,
            match_ic = 1.000,
            normalized_ic = 0.500,
            core_ic = 0.333,
            score = 0.500,
            query_information = 2.00,
            search_information = 4.00
        ),
        (
            query = raw"A..K",
            search = raw"^A.$",
            query_relationship = "Exact Overlap",
            search_relationship = "Exact Overlap",
            matched_pattern = raw"A.$",
            matched_positions = 1,
            match_ic = 1.000,
            normalized_ic = 0.500,
            core_ic = 0.500,
            score = 0.500,
            query_information = 2.00,
            search_information = 3.00
        ),
        (
            query = raw"[KR]L",
            search = raw"^AK$",
            query_relationship = "Degenerate Subsequence",
            search_relationship = "Variant Parent",
            matched_pattern = raw"[kr][$l]",
            matched_positions = 1,
            match_ic = 0.769,
            normalized_ic = 0.435,
            core_ic = 0.384,
            score = 0.435,
            query_information = 1.77,
            search_information = 4.00
        ),
        (
            query = raw"A[KR]L",
            search = raw"^AK$",
            query_relationship = "Degenerate Subsequence",
            search_relationship = "Variant Parent",
            matched_pattern = raw"A[kr][$l]",
            matched_positions = 2,
            match_ic = 1.769,
            normalized_ic = 0.639,
            core_ic = 0.590,
            score = 1.278,
            query_information = 2.77,
            search_information = 4.00
        ),
        (
            query = raw"A.K",
            search = raw"^AK$",
            query_relationship = "Degenerate Subsequence",
            search_relationship = "Variant Parent",
            matched_pattern = raw"[^a]aK",
            matched_positions = 1,
            match_ic = 1.000,
            normalized_ic = 0.500,
            core_ic = 0.333,
            score = 0.500,
            query_information = 2.00,
            search_information = 4.00
        ),
        (
            query = raw"A.K",
            search = raw"^A.$",
            query_relationship = "Exact Subsequence",
            search_relationship = "Exact Parent",
            matched_pattern = raw"A.[$k]",
            matched_positions = 1,
            match_ic = 1.000,
            normalized_ic = 0.500,
            core_ic = 0.500,
            score = 0.500,
            query_information = 2.00,
            search_information = 3.00
        )
    ]

    for case in oracle_cases
        result = compare(case.query, case.search, options)
        @test result.matched
        @test result.query_relationship == case.query_relationship
        @test result.search_relationship == case.search_relationship
        @test canonicalize_match_pattern(result.matched_pattern) ==
              canonicalize_match_pattern(case.matched_pattern)
        @test result.matched_positions == case.matched_positions
        @test result.match_ic ≈ case.match_ic atol = 1e-3
        @test result.normalized_ic ≈ case.normalized_ic atol = 1e-3
        @test result.core_ic ≈ case.core_ic atol = 1e-3
        @test result.score ≈ case.score atol = 1e-3
        @test result.query_information ≈ case.query_information atol = 1e-2
        @test result.search_information ≈ case.search_information atol = 1e-2
        @test !compare(case.query, case.search, default_options).matched
    end
end
