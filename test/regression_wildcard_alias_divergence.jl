@testitem "oracle wildcard-alias divergence fixture" begin
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

    function signature(result::ComparisonResult)
        return (
            matched = result.matched,
            normalized_query = result.normalized_query,
            normalized_search = result.normalized_search,
            query_relationship = result.query_relationship,
            search_relationship = result.search_relationship,
            matched_pattern = result.matched_pattern,
            matched_positions = result.matched_positions,
            match_ic = result.match_ic,
            normalized_ic = result.normalized_ic,
            score = result.score,
            query_information = result.query_information,
            search_information = result.search_information
        )
    end

    root = dirname(@__DIR__)
    motifs_path = joinpath(root, "data", "fixtures", "wildcard_alias_divergence_probe_motifs.tsv")
    oracle_path = joinpath(
        root, "data", "fixtures", "oracle_wildcard_alias_divergence_normalized.tsv")

    patterns = load_fixture_motifs(motifs_path)
    @test length(patterns) == 6

    oracle_rows = load_oracle_rows(oracle_path)
    observed_base_names = Set{String}()
    for row in oracle_rows
        push!(observed_base_names, base_name(row["Name1"]))
        push!(observed_base_names, base_name(row["Name2"]))
    end

    all_base_names = Set(["M0001", "M0002", "M0003", "M0004", "M0005", "M0006"])
    dropped_base_names = Set(["M0001", "M0003"])
    retained_base_names = Set(["M0002", "M0004", "M0005", "M0006"])
    @test issubset(dropped_base_names, setdiff(all_base_names, observed_base_names))
    @test issubset(retained_base_names, observed_base_names)

    options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)

    top_level_dot = compare("(Q|.)", "Q", options)
    top_level_x = compare("(Q|x)", "Q", options)
    top_level_X = compare("(Q|X)", "Q", options)
    top_level_reverse = compare("(.|Q)", "Q", options)
    @test signature(top_level_dot) == signature(top_level_x)
    @test signature(top_level_dot) == signature(top_level_X)
    @test signature(top_level_dot) == signature(top_level_reverse)
    @test top_level_dot.query_relationship == "Exact Match"
    @test top_level_dot.search_relationship == "Exact Match"
    @test top_level_dot.normalized_query == "(Q)|(.)"

    quantified_dot = compare("(Q|.){2}", "QQ", options)
    quantified_x = compare("(Q|x){2}", "QQ", options)
    quantified_X = compare("(Q|X){2}", "QQ", options)
    @test signature(quantified_dot) == signature(quantified_x)
    @test signature(quantified_dot) == signature(quantified_X)
    @test quantified_dot.query_relationship == "Exact Match"
    @test quantified_dot.search_relationship == "Exact Match"
    @test quantified_dot.normalized_query == "(Q{2})|(.{2})"

    embedded_dot = compare("A(Q|.)L", "AQL", options)
    embedded_x = compare("A(Q|x)L", "AQL", options)
    embedded_X = compare("A(Q|X)L", "AQL", options)
    @test signature(embedded_dot) == signature(embedded_x)
    @test signature(embedded_dot) == signature(embedded_X)
    @test embedded_dot.query_relationship == "Exact Match"
    @test embedded_dot.search_relationship == "Exact Match"
    @test embedded_dot.normalized_query == "(AQL)|(A.L)"

    embedded_dot_overlap = compare("A(Q|.)L", ".AA", options)
    embedded_x_overlap = compare("A(Q|x)L", ".AA", options)
    embedded_X_overlap = compare("A(Q|X)L", ".AA", options)
    @test signature(embedded_dot_overlap) == signature(embedded_x_overlap)
    @test signature(embedded_dot_overlap) == signature(embedded_X_overlap)
    @test embedded_dot_overlap.query_relationship == "Exact Overlap"
    @test embedded_dot_overlap.search_relationship == "Exact Overlap"
    @test embedded_dot_overlap.matched_pattern == "A"

    grouped_vs_explicit_overlap = compare("A(Q|.)L", ".[QL]", options)
    explicit_overlap = compare("(A.L)|(AQL)", ".[QL]", options)
    @test signature(grouped_vs_explicit_overlap) == signature(explicit_overlap)
    @test grouped_vs_explicit_overlap.normalized_query == "(AQL)|(A.L)"
    @test grouped_vs_explicit_overlap.query_relationship == "Variant Parent"
    @test grouped_vs_explicit_overlap.search_relationship == "Degenerate Subsequence"
    @test grouped_vs_explicit_overlap.matched_pattern == "q[ql]"

    grouped_specificity_overlap = compare("Q(.|[QL])", ".[QL]", options)
    explicit_specificity_overlap = compare("(Q[QL])|(Q.)", ".[QL]", options)
    @test signature(grouped_specificity_overlap) == signature(explicit_specificity_overlap)
    @test grouped_specificity_overlap.normalized_query == "(Q.)|(Q[QL])"
    @test grouped_specificity_overlap.query_relationship == "Variant Overlap"
    @test grouped_specificity_overlap.search_relationship == "Degenerate Overlap"
    @test grouped_specificity_overlap.matched_pattern == "[ql]"

    grouped_specificity_parent = compare("Q(.|[QL])", "H[LM]H(([KR][^H].)|(.[^H][KR]))", options)
    explicit_specificity_parent = compare("(Q[QL])|(Q.)", "H[LM]H(([KR][^H].)|(.[^H][KR]))", options)
    @test signature(grouped_specificity_parent) == signature(explicit_specificity_parent)
    @test grouped_specificity_parent.normalized_query == "(Q.)|(Q[QL])"
    @test grouped_specificity_parent.normalized_search ==
          "(H[LM]H.[ARNDCQEGILKMFPSTWYV][RK])|(H[LM]H[RK][ARNDCQEGILKMFPSTWYV].)"

    anchor_top_level = compare("(^Q)|(Q)", "Q", options)
    anchor_reverse = compare("(Q)|(^Q)", "Q", options)
    @test signature(anchor_top_level) == signature(anchor_reverse)
    @test anchor_top_level.normalized_query == "(Q)|(^Q)"
    @test anchor_top_level.query_relationship == "Exact Match"
    @test anchor_top_level.search_relationship == "Exact Match"
    @test anchor_top_level.matched_pattern == "Q"
end
