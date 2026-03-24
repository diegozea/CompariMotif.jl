@testitem "oracle corner-case parity fixture" begin
    function load_fixture_motifs(path::String)
        motifs = String[]
        for line in eachline(path)
            isempty(strip(line)) && continue
            startswith(line, '#') && continue
            cols = split(line, '\t')
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

        return (a["Name1"], a["Name2"], a["Motif1"], a["Motif2"]) <
               (b["Name1"], b["Name2"], b["Motif1"], b["Motif2"])
    end

    root = dirname(@__DIR__)
    motifs_path = joinpath(root, "data", "fixtures", "cornercase_probe_motifs.tsv")
    oracle_path = joinpath(root, "data", "fixtures", "oracle_cornercase_probe_normalized.tsv")

    raw_patterns = load_fixture_motifs(motifs_path)
    id_to_index = Dict("M$(lpad(string(idx), 4, '0'))" => idx
    for idx in eachindex(raw_patterns))

    oracle_rows = load_oracle_rows(oracle_path)
    observed_base_names = Set{String}()
    for row in oracle_rows
        push!(observed_base_names, base_name(row["Name1"]))
        push!(observed_base_names, base_name(row["Name2"]))
    end

    retained_indices = sort!(Int[id_to_index[base] for base in observed_base_names])
    retained_lookup = Dict(idx => retained_idx
    for (retained_idx, idx) in enumerate(retained_indices))
    retained_patterns = raw_patterns[retained_indices]

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    results = compare(retained_patterns, retained_patterns, options)

    best_by_pair = Dict{Tuple{Int, Int}, Dict{String, String}}()
    for row in oracle_rows
        qbase = base_name(row["Name1"])
        sbase = base_name(row["Name2"])
        qidx = id_to_index[qbase]
        sidx = id_to_index[sbase]
        qidx == sidx && continue
        pair = (retained_lookup[qidx], retained_lookup[sidx])
        if !haskey(best_by_pair, pair) || row_is_better(row, best_by_pair[pair])
            best_by_pair[pair] = row
        end
    end

    rejected_or_dropped = Set(["M$(lpad(string(idx), 4, '0'))"
                               for idx in eachindex(raw_patterns)
                               if !(idx in retained_indices)])
    @test rejected_or_dropped == Set(["M0007", "M0008", "M0014"])
    @test isempty(observed_base_names ∩ rejected_or_dropped)

    expected_pairs = Set(keys(best_by_pair))
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
