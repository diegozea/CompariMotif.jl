@testitem "oracle alternation variant-order parity fixture" begin
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

    function base_name(name::String)
        idx = findfirst(==('_'), name)
        return idx === nothing ? name : name[firstindex(name):prevind(name, idx)]
    end

    function row_is_better(a::Dict{String, String}, b::Dict{String, String})
        a_match_ic = parse(Float64, a["MatchIC"])
        b_match_ic = parse(Float64, b["MatchIC"])
        if a_match_ic != b_match_ic
            return a_match_ic > b_match_ic
        end

        a_match_pos = parse(Int, a["MatchPos"])
        b_match_pos = parse(Int, b["MatchPos"])
        if a_match_pos != b_match_pos
            return a_match_pos > b_match_pos
        end

        a_score = parse(Float64, a["Score"])
        b_score = parse(Float64, b["Score"])
        if a_score != b_score
            return a_score > b_score
        end

        return (a["Name1"], a["Name2"], a["Motif1"], a["Motif2"], a["Match"]) <
               (b["Name1"], b["Name2"], b["Motif1"], b["Motif2"], b["Match"])
    end

    root = dirname(@__DIR__)
    motifs_path = joinpath(root, "data", "fixtures", "variant_order_probe_motifs.tsv")
    oracle_path = joinpath(root, "data", "fixtures", "oracle_variant_order_probe_normalized.tsv")

    patterns = load_fixture_motifs(motifs_path)
    id_to_index = Dict("M$(lpad(string(idx), 4, '0'))" => idx
    for idx in eachindex(patterns))
    @test length(id_to_index) == length(patterns)

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    results = compare(patterns, patterns, options)
    oracle_rows = load_oracle_rows(oracle_path)

    best_by_pair = Dict{Tuple{Int, Int}, Dict{String, String}}()
    for row in oracle_rows
        qi = id_to_index[base_name(row["Name1"])]
        si = id_to_index[base_name(row["Name2"])]
        qi == si && continue
        pair = (qi, si)
        if !haskey(best_by_pair, pair) || row_is_better(row, best_by_pair[pair])
            best_by_pair[pair] = row
        end
    end

    for (pair, row) in best_by_pair
        qi, si = pair
        result = results[qi, si]
        @test result.matched
        @test result.query_relationship == row["Sim1"]
        @test result.search_relationship == row["Sim2"]
        @test result.matched_pattern == row["Match"]
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
    @test observed_pairs == Set(keys(best_by_pair))

    @test results[1, 2].query_relationship == "Exact Overlap"
    @test results[1, 2].search_relationship == "Exact Overlap"
    @test results[1, 2].matched_pattern == "QQ"

    @test results[1, 3].query_relationship == "Variant Match"
    @test results[1, 3].search_relationship == "Degenerate Match"
    @test results[1, 3].matched_pattern == "[nq]Q[st]"

    @test results[2, 1].query_relationship == "Exact Overlap"
    @test results[2, 1].search_relationship == "Exact Overlap"
    @test results[3, 1].query_relationship == "Degenerate Match"
    @test results[3, 1].search_relationship == "Variant Match"

    @test results[4, 5].query_relationship == "Exact Match"
    @test results[4, 5].search_relationship == "Exact Match"
    @test results[4, 5].matched_pattern == "QQ"

    @test results[5, 4].query_relationship == "Exact Match"
    @test results[5, 4].search_relationship == "Exact Match"
    @test results[5, 4].matched_pattern == "QQ"

    @test results[6, 7].query_relationship == "Exact Parent"
    @test results[6, 7].search_relationship == "Exact Subsequence"
    @test results[6, 7].matched_pattern == "AA"

    @test results[7, 6].query_relationship == "Exact Subsequence"
    @test results[7, 6].search_relationship == "Exact Parent"
    @test results[7, 6].matched_pattern == "AA"

    @test results[8, 9].query_relationship == "Variant Match"
    @test results[8, 9].search_relationship == "Degenerate Match"
    @test results[8, 9].matched_pattern == "[aq]"

    @test results[10, 9].query_relationship == "Variant Match"
    @test results[10, 9].search_relationship == "Degenerate Match"
    @test results[10, 9].matched_pattern == "[aq]"
end
