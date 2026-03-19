@testitem "oracle matchfix parity fixtures" begin
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
                lowercase_body = all(islowercase, body)
                normalized_body = uppercase(String(body))
                mask = CompariMotif._class_mask(normalized_body, spec)
                print(
                    io,
                    CompariMotif._mask_to_symbol(
                        mask,
                        spec;
                        as_lowercase = lowercase_body,
                        wildcard_symbol = "x"
                    )
                )
                i = nextind(pattern, close_idx)
                continue
            end
            print(io, pattern[i])
            i = nextind(pattern, i)
        end
        return String(take!(io))
    end

    root = dirname(@__DIR__)
    motifs_path = joinpath(root, "data", "fixtures", "matchfix_probe_motifs.tsv")
    patterns = load_fixture_motifs(motifs_path)
    motif_to_index = Dict(motif => idx for (idx, motif) in enumerate(patterns))
    @test length(motif_to_index) == length(patterns)

    fixture_specs = [
        (
            oracle_path = joinpath(
                root, "data", "fixtures", "oracle_query_fixed_probe_normalized.tsv"),
            options = ComparisonOptions(;
                min_shared_positions = 1,
                normalized_ic_cutoff = 0.0,
                matchfix = :query_fixed
            ),
            absent_pairs = [
                ("[IPV]IAL", "IPI"),
                ("A[KR]L", "RL"),
                ("K.[KR].", ".[KR]"),
                ("A.K", ".K"),
                ("A.K", "A."),
                ("AK.A", "AK.")
            ],
            present_pairs = [
                ("AKLI", "KLI"),
                ("[IPV]IAL", "IAL"),
                ("A[KR]L", "[KR]L"),
                (".[KR]", "K.[KR]."),
                (".AA", ".A"),
                ("AA.", "A."),
                ("AK.A", "K.A"),
                ("A.KA", "A.K"),
                (".K", "A.K")
            ]
        ),
        (
            oracle_path = joinpath(
                root, "data", "fixtures", "oracle_search_fixed_probe_normalized.tsv"),
            options = ComparisonOptions(;
                min_shared_positions = 1,
                normalized_ic_cutoff = 0.0,
                matchfix = :search_fixed
            ),
            absent_pairs = [
                ("IPI", "[IPV]IAL"),
                ("RL", "A[KR]L"),
                (".[KR]", "K.[KR]."),
                (".K", "A.K"),
                ("A.", "A.K"),
                ("AK.", "AK.A")
            ],
            present_pairs = [
                ("KLI", "AKLI"),
                ("IAL", "[IPV]IAL"),
                ("[KR]L", "A[KR]L"),
                ("K.[KR].", ".[KR]"),
                (".A", ".AA"),
                ("A.", "AA."),
                ("K.A", "AK.A"),
                ("A.K", "A.KA"),
                ("A.K", ".K")
            ]
        ),
        (
            oracle_path = joinpath(
                root, "data", "fixtures", "oracle_both_fixed_probe_normalized.tsv"),
            options = ComparisonOptions(;
                min_shared_positions = 1,
                normalized_ic_cutoff = 0.0,
                matchfix = :both_fixed
            ),
            absent_pairs = [
                ("AKLI", "KLIA"),
                ("[IPV]IAL", "IPI"),
                ("A[KR]L", "RL"),
                ("RL", "A[KR]L"),
                (".[KR]", "K.[KR]."),
                ("K.[KR].", ".[KR]"),
                ("A.K", ".K"),
                (".K", "A.K"),
                ("A.K", "A."),
                ("A.", "A.K"),
                ("AK.A", "AK."),
                ("AK.", "AK.A")
            ],
            present_pairs = [
                ("AKLI", "KLI"),
                ("[IPV]IAL", "IAL"),
                ("A[KR]L", "[KR]L"),
                ("[KR]L", "A[KR]L"),
                (".AA", ".A"),
                (".A", ".AA"),
                ("AA.", "A."),
                ("A.", "AA."),
                ("AK.A", "K.A"),
                ("K.A", "AK.A"),
                ("A.KA", "A.K"),
                ("A.K", "A.KA")
            ]
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
            @test result.query_relationship == row["Sim1"]
            @test result.search_relationship == row["Sim2"]
            @test canonicalize_match_pattern(result.matched_pattern) ==
                  canonicalize_match_pattern(replace(row["Match"], '.' => 'x'))
            @test result.matched_positions == parse(Int, row["MatchPos"])
            @test result.match_ic ≈ parse(Float64, row["MatchIC"]) atol = 1e-3
            @test result.normalized_ic ≈ parse(Float64, row["NormIC"]) atol = 1e-3
            @test result.core_ic ≈ parse(Float64, row["CoreIC"]) atol = 1e-3
            @test result.score ≈ parse(Float64, row["Score"]) atol = 1e-3
            @test result.query_information ≈ parse(Float64, row["Info1"]) atol = 1e-2
            @test result.search_information ≈ parse(Float64, row["Info2"]) atol = 1e-2
        end

        observed_pairs = Set{Tuple{Int, Int}}()
        for i in axes(results, 1), j in axes(results, 2)

            if i != j && results[i, j].matched
                push!(observed_pairs, (i, j))
            end
        end
        @test observed_pairs == expected_pairs

        for (query, search) in spec.absent_pairs
            @test !results[motif_to_index[query], motif_to_index[search]].matched
        end
        for (query, search) in spec.present_pairs
            @test results[motif_to_index[query], motif_to_index[search]].matched
        end
    end
end
