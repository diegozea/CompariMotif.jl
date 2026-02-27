using TestItems

@testitem "oracle defaults parity fixture" begin
    using Test
    using CompariMotif

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

    oracle_to_package(rel::String) = replace(rel, "Ugly" => "Complex")

    root = dirname(@__DIR__)
    motifs_path = joinpath(root, "data", "fixtures", "default_probe_motifs.tsv")
    oracle_path = joinpath(root, "data", "fixtures", "oracle_default_probe_normalized.tsv")

    patterns = load_fixture_motifs(motifs_path)
    motif_to_index = Dict(motif => idx for (idx, motif) in enumerate(patterns))
    @test length(motif_to_index) == length(patterns)

    results = compare(patterns, patterns, ComparisonOptions())
    oracle_rows = load_oracle_rows(oracle_path)

    expected_pairs = Set{Tuple{Int, Int}}()
    for row in oracle_rows
        qi = motif_to_index[row["Motif1"]]
        si = motif_to_index[row["Motif2"]]
        push!(expected_pairs, (qi, si))

        result = results[qi, si]
        is_ugly = occursin("Ugly", row["Sim1"]) || occursin("Ugly", row["Sim2"])

        @test result.matched
        @test result.query_relationship == oracle_to_package(row["Sim1"])
        @test result.search_relationship == oracle_to_package(row["Sim2"])
        @test result.matched_positions == parse(Int, row["MatchPos"])
        if !is_ugly
            @test result.match_ic ≈ parse(Float64, row["MatchIC"]) atol = 1e-3
            @test result.normalized_ic ≈ parse(Float64, row["NormIC"]) atol = 1e-3
            @test result.score ≈ parse(Float64, row["Score"]) atol = 1e-3
            @test result.query_information ≈ parse(Float64, row["Info1"]) atol = 1e-2
            @test result.search_information ≈ parse(Float64, row["Info2"]) atol = 1e-2
        end
    end

    observed_pairs = Set{Tuple{Int, Int}}()
    for i in axes(results, 1), j in axes(results, 2)

        if i != j && results[i, j].matched
            push!(observed_pairs, (i, j))
        end
    end
    @test observed_pairs == expected_pairs
end
