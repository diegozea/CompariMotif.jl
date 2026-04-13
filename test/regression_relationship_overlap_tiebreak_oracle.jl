@testitem "oracle relationship-overlap tie-break parity fixture" begin
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

    function canonicalize_match_pattern(pattern::AbstractString)
        spec = CompariMotif._alphabet_spec(ProteinAlphabet())
        io = IOBuffer()
        i = firstindex(pattern)
        while i <= lastindex(pattern)
            if pattern[i] == '['
                close_idx = findnext(==(']'), pattern, nextind(pattern, i))
                close_idx === nothing && error("Malformed match pattern: $pattern")
                body = pattern[nextind(pattern, i):prevind(pattern, close_idx)]
                negated = !isempty(body) &&
                          first(body) == '^' &&
                          any(c -> haskey(spec.index, uppercase(c)),
                              body[nextind(body, firstindex(body)):end])
                anchors = IOBuffer()
                residue_chars = Char[]
                for idx in eachindex(body)
                    char = body[idx]
                    if negated && idx == firstindex(body)
                        continue
                    elseif char == '^' || char == '$'
                        print(anchors, char)
                    else
                        push!(residue_chars, char)
                    end
                end
                anchor_text = String(take!(anchors))
                residues = String(residue_chars)
                normalized_residues = if isempty(residues)
                    ""
                else
                    mask = if negated
                        CompariMotif._mask_complement(
                            CompariMotif._class_mask(uppercase(residues), spec),
                            spec
                        )
                    else
                        CompariMotif._class_mask(uppercase(residues), spec)
                    end
                    symbol = CompariMotif._mask_to_symbol(
                        mask,
                        spec;
                        as_lowercase = all(islowercase, residues)
                    )
                    startswith(symbol, '[') && endswith(symbol, ']') ? symbol[2:(end - 1)] :
                    symbol
                end
                print(io, '[', anchor_text, normalized_residues, ']')
                i = nextind(pattern, close_idx)
                continue
            end
            print(io, pattern[i])
            i = nextind(pattern, i)
        end
        return String(take!(io))
    end

    oracle_to_package(rel::String) = replace(rel, "Ugly" => "Complex")

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

        return (a["Name1"], a["Name2"], a["Motif1"],
            a["Motif2"], a["Sim1"], a["Sim2"], a["Match"]) <
               (b["Name1"], b["Name2"], b["Motif1"],
            b["Motif2"], b["Sim1"], b["Sim2"], b["Match"])
    end

    root = dirname(@__DIR__)
    motifs_path = joinpath(
        root,
        "data",
        "fixtures",
        "relationship_overlap_tiebreak_probe_motifs.tsv"
    )
    oracle_path = joinpath(
        root,
        "data",
        "fixtures",
        "oracle_relationship_overlap_tiebreak_probe_normalized.tsv"
    )

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

    target_pairs = Set([
        (1, 2), (2, 1),
        (3, 4), (4, 3),
        (5, 6), (6, 5),
        (7, 8), (8, 7),
        (9, 10), (10, 9),
        (11, 12),
        (2, 13), (13, 2),
        (4, 14), (14, 4),
        (5, 8), (8, 5),
        (1, 15), (15, 1),
        (14, 16), (16, 14)
    ])
    for (pair, row) in best_by_pair
        pair in target_pairs || continue
        qi, si = pair
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
    end

    @test results[1, 2].query_relationship == "Complex Parent"
    @test results[1, 2].search_relationship == "Complex Subsequence"
    @test results[1, 2].matched_pattern == "[as].s\$"

    @test results[3, 4].query_relationship == "Complex Subsequence"
    @test results[3, 4].search_relationship == "Complex Parent"

    @test results[5, 6].query_relationship == "Complex Subsequence"
    @test results[5, 6].search_relationship == "Complex Parent"
    @test results[5, 6].matched_pattern == "v.g[RK][rk]"

    @test results[7, 8].query_relationship == "Complex Match"
    @test results[7, 8].search_relationship == "Complex Match"

    @test results[9, 10].query_relationship == "Complex Subsequence"
    @test results[9, 10].search_relationship == "Complex Parent"
    @test results[10, 9].query_relationship == "Complex Parent"
    @test results[10, 9].search_relationship == "Complex Subsequence"

    @test results[11, 12].query_relationship == "Complex Subsequence"
    @test results[11, 12].search_relationship == "Complex Parent"

    @test results[13, 2].query_relationship == "Complex Overlap"
    @test results[13, 2].search_relationship == "Complex Overlap"
    @test results[13, 2].matched_pattern == "[as]s"
    @test results[2, 13].query_relationship == "Complex Overlap"
    @test results[2, 13].search_relationship == "Complex Overlap"

    @test results[14, 4].query_relationship == "Complex Overlap"
    @test results[14, 4].search_relationship == "Complex Overlap"
    @test results[14, 4].matched_pattern == "h[lm]"
    @test results[4, 14].query_relationship == "Degenerate Overlap"
    @test results[4, 14].search_relationship == "Variant Overlap"
    @test results[4, 14].matched_pattern == "[rk]"

    @test results[5, 8].query_relationship == "Complex Subsequence"
    @test results[5, 8].search_relationship == "Complex Parent"
    @test results[8, 5].query_relationship == "Complex Parent"
    @test results[8, 5].search_relationship == "Complex Subsequence"

    @test results[1, 15].query_relationship == "Complex Parent"
    @test results[1, 15].search_relationship == "Complex Subsequence"
    @test results[1, 15].matched_pattern == "Aaka"
    @test results[15, 1].query_relationship == "Exact Overlap"
    @test results[15, 1].search_relationship == "Exact Overlap"
    @test results[15, 1].matched_pattern == "A"

    @test results[14, 16].query_relationship == "Complex Subsequence"
    @test results[14, 16].search_relationship == "Complex Parent"
    @test results[14, 16].matched_pattern == "[rk][rk]l"
    @test results[16, 14].query_relationship == "Complex Parent"
    @test results[16, 14].search_relationship == "Complex Subsequence"
    @test results[16, 14].matched_pattern == "[rk]wl"
end
