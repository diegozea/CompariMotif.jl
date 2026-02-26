"""
    normalize_motif(motif::AbstractString; alphabet::Symbol = :protein) -> String

Parse and canonicalize a motif expression into a deterministic representation.
Supported syntax includes fixed residues, bracket classes (including negation),
`x`/`.` wildcards, `^`/`\$` termini, and `{n}`/`{m,n}` (or `(n)`/`(m,n)`) repeat quantifiers.

$(_DOC_OPTIONS_REF)
$(_DOC_COMPARE_REF)
"""
function normalize_motif(motif::AbstractString; alphabet::Symbol = :protein)
    # Reuse the normal parser pipeline with permissive thresholds.
    # `normalize_motif` only needs deterministic parsing/canonicalization.
    options = ComparisonOptions(; alphabet, min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    return _parse_motif(motif, options).normalized
end

"""
    compare(a::AbstractString, b::AbstractString, options::ComparisonOptions) -> ComparisonResult
    compare(motifs::AbstractVector{<:AbstractString},
            db::AbstractVector{<:AbstractString},
            options::ComparisonOptions) -> Matrix{ComparisonResult}
    compare(motifs::AbstractVector{<:AbstractString},
            options::ComparisonOptions) -> Matrix{ComparisonResult}

Compare motifs according to the CompariMotif scoring scheme described in
Edwards et al. (2008).

- Pairwise mode compares one query motif against one search motif.
- Matrix mode computes all pairwise query-vs-database comparisons.
- All-vs-all mode is a convenience alias for `compare(motifs, motifs, options)`.

$(_DOC_OPTIONS_REF)
$(_DOC_VARIANT_SIZE_REF)
$(_DOC_RESULT_REF)
$(_DOC_NORMALIZE_REF)
$(_DOC_TSV_REF)
"""
function compare end

function compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)
    # Parse both motifs once, then run the shared comparison core.
    parsed_a = _parse_motif(a, options)
    parsed_b = _parse_motif(b, options)
    return _compare_parsed(parsed_a, parsed_b, options)
end

function compare(
        motifs::AbstractVector{<:AbstractString},
        db::AbstractVector{<:AbstractString},
        options::ComparisonOptions
)
    # Parse input vectors upfront so malformed motifs fail early.
    parsed_queries = [_parse_motif(motif, options) for motif in motifs]
    parsed_db = [_parse_motif(motif, options) for motif in db]

    results = Matrix{ComparisonResult}(undef, length(parsed_queries), length(parsed_db))
    # Dense pairwise evaluation in deterministic row/column order.
    for i in eachindex(parsed_queries)
        for j in eachindex(parsed_db)
            results[i, j] = _compare_parsed(parsed_queries[i], parsed_db[j], options)
        end
    end
    return results
end

function compare(motifs::AbstractVector{<:AbstractString}, options::ComparisonOptions)
    # Convenience all-vs-all matrix mode.
    compare(motifs, motifs, options)
end
