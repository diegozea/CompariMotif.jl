"""
    write_results_tsv(path, motifs, db, results)

Write pairwise comparison results to a deterministic TSV file.
Pass `results` from `compare(motifs, db, options)` with matching matrix dimensions.

$(_DOC_COMPARE_REF)
$(_DOC_RESULT_REF)
"""
function write_results_tsv(path::AbstractString,
        motifs::AbstractVector{<:AbstractString},
        db::AbstractVector{<:AbstractString},
        results::Matrix{ComparisonResult})
    # Guard against mismatched caller inputs before writing any file content.
    size(results) == (length(motifs), length(db)) ||
        throw(ArgumentError("`results` dimensions must match `(length(motifs), length(db))`."))

    open(path, "w") do io
        # Stable column layout used by fixtures and downstream scripts.
        println(io,
            join(
                [
                    "QueryIndex",
                    "SearchIndex",
                    "QueryMotif",
                    "SearchMotif",
                    "Matched",
                    "QueryRelationship",
                    "SearchRelationship",
                    "QueryCode",
                    "SearchCode",
                    "MatchedPattern",
                    "MatchPos",
                    "MatchIC",
                    "NormIC",
                    "CoreIC",
                    "Score",
                    "InfoQuery",
                    "InfoSearch"
                ],
                '\t'))
        for i in eachindex(motifs)
            for j in eachindex(db)
                result = results[i, j]
                # Emit one deterministic row per matrix cell.
                println(io,
                    join(
                        [
                            string(i),
                            string(j),
                            motifs[i],
                            db[j],
                            result.matched ? "1" : "0",
                            result.query_relationship,
                            result.search_relationship,
                            result.query_relationship_code,
                            result.search_relationship_code,
                            result.matched_pattern,
                            string(result.matched_positions),
                            # Keep numeric formatting fixed to aid regression diffs.
                            string(round(result.match_ic; digits = 6)),
                            string(round(result.normalized_ic; digits = 6)),
                            string(round(result.core_ic; digits = 6)),
                            string(round(result.score; digits = 6)),
                            string(round(result.query_information; digits = 6)),
                            string(round(result.search_information; digits = 6))
                        ],
                        '\t'))
            end
        end
    end
    return path
end
