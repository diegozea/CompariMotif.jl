"""
    _normalize_motif(motif::AbstractString; alphabet = ProteinAlphabet())::String

Internal parser helper that renders a motif into the deterministic canonical
syntax used by the comparison pipeline.
"""
function _normalize_motif(motif::AbstractString; alphabet::_AlphabetValue = ProteinAlphabet())
    # Reuse the normal parser pipeline with permissive thresholds.
    # `_normalize_motif` only needs deterministic parsing/canonicalization.
    options = ComparisonOptions(; alphabet, min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    return _parse_motif(motif, options).normalized
end

"""
    compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)::ComparisonResult
    compare(motif::AbstractString,
            db::AbstractVector{<:AbstractString},
            options::ComparisonOptions)::Matrix{ComparisonResult}
    compare(motifs::AbstractVector{<:AbstractString},
            db::AbstractVector{<:AbstractString},
            options::ComparisonOptions)::Matrix{ComparisonResult}
    compare(motifs::AbstractVector{<:AbstractString},
            options::ComparisonOptions)::Matrix{ComparisonResult}

Compare motifs according to the CompariMotif scoring scheme described by
[Edwards2008CompariMotif](@citet).

- Pairwise mode compares one query motif against one target motif.
- Query-search mode compares one query motif against a database of target motifs.
- Database mode compares one motif collection against another.
- Passing a single motif collection allows for all-vs-all comparison within that collection.

# Examples
```jldoctest
julia> using CompariMotif

julia> options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> result = compare("RKLI", "R[KR]L[IV]", options);

julia> result.matched
true
```

# References
* [Edwards2008CompariMotif](@cite) Edwards et al. Bioinformatics 24(10):1307-1309 (2008)

$(_DOC_OPTIONS_REF)
$(_DOC_VARIANT_SIZE_REF)
$(_DOC_RESULT_REF)
$(_DOC_TABLE_REF)
"""
function compare end

"""
    compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)::ComparisonResult

Compare one query motif against one target motif.
"""
function compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)
    # Parse both motifs once, then run the shared comparison core.
    parsed_a = _parse_motif(a, options)
    parsed_b = _parse_motif(b, options)
    return _compare_parsed(parsed_a, parsed_b, options)
end

"""
    compare(motif, db, options)::Matrix{ComparisonResult}

Search one query motif against a database of target motifs.
"""
function compare(
        motif::AbstractString,
        db::AbstractVector{<:AbstractString},
        options::ComparisonOptions
)
    compare([motif], db, options)
end

"""
    compare(motifs, db, options)::Matrix{ComparisonResult}

Compare one motif collection against another motif collection.
"""
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

"""
    compare(motifs, options)::Matrix{ComparisonResult}

Convenience all-vs-all comparison of one motif collection.
"""
function compare(motifs::AbstractVector{<:AbstractString}, options::ComparisonOptions)
    # Convenience all-vs-all comparison mode.
    compare(motifs, motifs, options)
end
