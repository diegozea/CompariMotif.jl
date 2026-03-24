@testitem "oracle DNA fixtures" begin
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

    function load_residue_frequencies(path::String)
        freqs = Dict{Char, Float64}()
        for line in eachline(path)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, '#') && continue
            stripped == "AA\tFREQ" && continue
            cols = split(stripped, '\t')
            length(cols) == 2 || error("Malformed frequency row: $line")
            cols[1] == "Total" && continue
            length(cols[1]) == 1 || error("Invalid residue key in row: $line")
            freqs[only(cols[1])] = parse(Float64, cols[2])
        end
        return freqs
    end

    oracle_to_package(rel::String) = replace(rel, "Ugly" => "Complex")

    function assert_fixture(motifs_path::String, oracle_path::String, options::ComparisonOptions)
        patterns = load_fixture_motifs(motifs_path)
        motif_to_index = Dict(motif => idx for (idx, motif) in enumerate(patterns))
        @test length(motif_to_index) == length(patterns)

        results = compare(patterns, patterns, options)
        oracle_rows = load_oracle_rows(oracle_path)
        oracle_pairs = Set((row["Motif1"], row["Motif2"]) for row in oracle_rows)

        # `ATG` and `[AGT]TG` are the published DNA example pair used in the
        # README and external API docs; the remaining motifs are fixture-only
        # corner cases for wildcard, negation, overlap, and anchor behavior.
        documented_pairs = Set([("ATG", "[AGT]TG"), ("[AGT]TG", "ATG")])
        @test documented_pairs ⊆ oracle_pairs

        expected_pairs = Set{Tuple{Int, Int}}()
        for row in oracle_rows
            qi = motif_to_index[row["Motif1"]]
            si = motif_to_index[row["Motif2"]]
            push!(expected_pairs, (qi, si))

            result = results[qi, si]
            @test result.matched
            @test result.query_relationship == oracle_to_package(row["Sim1"])
            @test result.search_relationship == oracle_to_package(row["Sim2"])
            @test result.matched_pattern == row["Match"]
            @test result.matched_positions == parse(Int, row["MatchPos"])
            @test result.match_ic ≈ parse(Float64, row["MatchIC"]) atol = 1e-3
            @test result.normalized_ic ≈ parse(Float64, row["NormIC"]) atol = 1e-3
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
    end

    root = dirname(@__DIR__)
    motifs_path = joinpath(root, "data", "fixtures", "dna_probe_motifs.tsv")
    uniform_oracle_path = joinpath(root, "data", "fixtures", "oracle_dna_probe_normalized.tsv")
    nonuniform_oracle_path = joinpath(
        root, "data", "fixtures", "oracle_dna_nonuniform_probe_normalized.tsv"
    )
    frequency_path = joinpath(root, "data", "fixtures", "dna_nonuniform_probe.aafreq.tsv")

    @testset "uniform DNA frequencies" begin
        options = ComparisonOptions(;
            alphabet = DNAAlphabet(),
            min_shared_positions = 1,
            normalized_ic_cutoff = 0.0
        )
        assert_fixture(motifs_path, uniform_oracle_path, options)
    end

    @testset "non-uniform DNA frequencies" begin
        options = ComparisonOptions(;
            alphabet = DNAAlphabet(),
            residue_frequencies = load_residue_frequencies(frequency_path),
            min_shared_positions = 1,
            normalized_ic_cutoff = 0.0
        )
        assert_fixture(motifs_path, nonuniform_oracle_path, options)
    end
end
