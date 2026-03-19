"""
    _match_symbol(qpos, spos, intersection, relation, mismatch, spec)::String

Render one output symbol for the overlap pattern.
"""
function _match_symbol(
        qpos::_Position,
        spos::_Position,
        intersection::ResidueMask,
        relation::_RelationshipType,
        mismatch::Bool,
        spec::_AlphabetSpec
)
    if qpos.kind == _NTERMINUS
        return "^"
    elseif qpos.kind == _CTERMINUS
        return "\$"
    end

    qwild = _is_wildcard(qpos, spec.mask)
    swild = _is_wildcard(spos, spec.mask)
    if qwild && swild
        # Preserve compact wildcard representation in overlap output.
        return "x"
    end

    qclass = ResidueClass(qpos.mask)
    sclass = ResidueClass(spos.mask)
    if mismatch
        # Mismatch still contributes a descriptive symbol in the matched pattern.
        union_mask = unionclass(qclass, sclass).mask
        return _mask_to_symbol(union_mask, spec; as_lowercase = true, wildcard_symbol = "x")
    end
    if relation == _REL_EXACT
        return _mask_to_symbol(intersection, spec; as_lowercase = false, wildcard_symbol = "x")
    end
    if qwild
        return _mask_to_symbol(spos.mask, spec; as_lowercase = true, wildcard_symbol = "x")
    elseif swild
        return _mask_to_symbol(qpos.mask, spec; as_lowercase = true, wildcard_symbol = "x")
    end
    union_mask = unionclass(qclass, sclass).mask
    return _mask_to_symbol(union_mask, spec; as_lowercase = true, wildcard_symbol = "x")
end

"""
    _compare_positions(qpos, spos, options, spec)

Compare one query/search position pair and return matching diagnostics.
"""
function _compare_positions(
        qpos::_Position,
        spos::_Position,
        options::ComparisonOptions,
        spec::_AlphabetSpec
)
    return _compare_positions(
        qpos, spos, options, spec, _residue_frequency_vector(options, spec))
end

function _compare_positions(
        qpos::_Position,
        spos::_Position,
        options::ComparisonOptions,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    # Anchors are valid only against the same anchor type.
    if qpos.kind !== _RESIDUE || spos.kind !== _RESIDUE
        if qpos.kind == spos.kind
            return (hard_mismatch = false, mismatch = false,
                relation = _REL_EXACT, intersection = ResidueMask(0), ic = 1.0,
                core_ic_denominator = 1.0,
                contributes_position = true, exact_fixed = false)
        end
        return (hard_mismatch = true, mismatch = true,
            relation = _REL_COMPLEX, intersection = ResidueMask(0), ic = 0.0,
            core_ic_denominator = 0.0,
            contributes_position = false, exact_fixed = false)
    end

    q_ic = _position_ic(qpos, spec, residue_frequencies)
    s_ic = _position_ic(spos, spec, residue_frequencies)

    qclass = ResidueClass(qpos.mask)
    sclass = ResidueClass(spos.mask)
    # Residue match operates on mask set intersection.
    intersection = qclass.mask & sclass.mask
    if !overlaps(qclass, sclass)
        # No shared residue -> mismatch position (allowed only if mismatch budget allows).
        return (hard_mismatch = false, mismatch = true,
            relation = _REL_COMPLEX, intersection = intersection, ic = 0.0,
            core_ic_denominator = max(q_ic, s_ic),
            contributes_position = false, exact_fixed = false)
    end

    # Determine relationship by set equality/subset/superset/partial overlap.
    relation = if qclass.mask == sclass.mask
        _REL_EXACT
    elseif is_subset(qclass, sclass)
        _REL_VARIANT
    elseif is_subset(sclass, qclass)
        _REL_DEGENERATE
    elseif options.allow_ambiguous_overlap
        _REL_COMPLEX
    else
        return (hard_mismatch = true, mismatch = true, relation = _REL_COMPLEX,
            intersection = intersection, ic = 0.0,
            core_ic_denominator = max(q_ic, s_ic),
            contributes_position = false, exact_fixed = false)
    end

    # Partially overlapping ambiguous classes use the union class IC, matching
    # the oracle's per-position "Ugly" scoring. Subset/superset matches still
    # use the less-specific side.
    ic = if relation == _REL_COMPLEX
        _position_ic(_Position(_RESIDUE, unionclass(qclass, sclass).mask), spec, residue_frequencies)
    else
        min(q_ic, s_ic)
    end
    # Matched position count excludes wildcard-vs-wildcard.
    contributes = !_is_wildcard(qpos, spec.mask) && !_is_wildcard(spos, spec.mask)
    core_ic_denominator = max(q_ic, s_ic)
    # Used for matchfix tie-break logic.
    exact_fixed = relation == _REL_EXACT && _is_fixed(qpos) && _is_fixed(spos)
    return (
        hard_mismatch = false,
        mismatch = false,
        relation = relation,
        intersection = intersection,
        ic = ic,
        core_ic_denominator = core_ic_denominator,
        contributes_position = contributes,
        exact_fixed = exact_fixed
    )
end

"""
    _query_fixed_required(mode::Symbol)::Bool

Return `true` when query fixed residues must match exactly.
"""
function _query_fixed_required(mode::Symbol)
    # Query fixed constraints apply in QueryFixed and BothFixed modes.
    if mode === :query_fixed || mode === :both_fixed
        return true
    elseif mode === :none || mode === :search_fixed
        return false
    end
    throw(_matchfix_argument_error())
end

"""
    _search_fixed_required(mode::Symbol)::Bool

Return `true` when search fixed residues must match exactly.
"""
function _search_fixed_required(mode::Symbol)
    # Search fixed constraints apply in SearchFixed and BothFixed modes.
    if mode === :search_fixed || mode === :both_fixed
        return true
    elseif mode === :none || mode === :query_fixed
        return false
    end
    throw(_matchfix_argument_error())
end

"""
    _clips_fixed_prefix(positions, overlap_start)::Bool

Return `true` when the constrained side clips one or more fixed residues before
the aligned span.
"""
function _clips_fixed_prefix(
        positions::AbstractVector{_Position},
        overlap_start::Int
)
    @inbounds for idx in firstindex(positions):(overlap_start - 1)
        if _is_fixed(positions[idx])
            return true
        end
    end
    return false
end

"""
    _clips_fixed_suffix(positions, overlap_end)::Bool

Return `true` when the constrained side clips one or more fixed residues after
the aligned span.
"""
function _clips_fixed_suffix(
        positions::AbstractVector{_Position},
        overlap_end::Int
)
    @inbounds for idx in (overlap_end + 1):lastindex(positions)
        if _is_fixed(positions[idx])
            return true
        end
    end
    return false
end

"""
    _is_exact_contained_alignment(query_positions, search_positions, overlap_start, search_overlap_start, overlap_length)::Bool

Return `true` when the aligned span covers the shorter motif completely and the
contained alignment is exact position-for-position. The oracle only exempts
these exact subsequence/parent cases from clipped fixed-position rejection
under `matchfix` when the shorter motif keeps an informative boundary residue on
each side where the constrained longer motif clips fixed positions.
"""
function _is_exact_contained_alignment(
        query_positions::AbstractVector{_Position},
        search_positions::AbstractVector{_Position},
        overlap_start::Int,
        search_overlap_start::Int,
        overlap_length::Int,
        query_clips_fixed_prefix::Bool,
        query_clips_fixed_suffix::Bool,
        search_clips_fixed_prefix::Bool,
        search_clips_fixed_suffix::Bool,
        wildcard_mask::ResidueMask
)
    qlen = length(query_positions)
    slen = length(search_positions)
    overlap_length == min(qlen, slen) || return false

    if qlen < slen
        return _matches_matchfix_exact_subsequence(
            query_positions,
            search_positions,
            search_overlap_start,
            wildcard_mask,
            search_clips_fixed_prefix,
            search_clips_fixed_suffix
        )
    elseif slen < qlen
        return _matches_matchfix_exact_subsequence(
            search_positions,
            query_positions,
            overlap_start,
            wildcard_mask,
            query_clips_fixed_prefix,
            query_clips_fixed_suffix
        )
    end
    return _matches_exact_subsequence(query_positions, search_positions, 1)
end

"""
    _evaluate_alignment(query_variant, search_variant, shift, options, spec)

Evaluate one concrete shift between two expanded motif variants.
Returns `_Candidate` when all thresholds pass, otherwise `nothing`.
"""
function _evaluate_alignment(
        query_variant::_MotifVariant,
        search_variant::_MotifVariant,
        shift::Int,
        options::ComparisonOptions,
        spec::_AlphabetSpec
)
    return _evaluate_alignment(
        query_variant,
        search_variant,
        shift,
        options,
        spec,
        _residue_frequency_vector(options, spec)
    )
end

function _evaluate_alignment(
        query_variant::_MotifVariant,
        search_variant::_MotifVariant,
        shift::Int,
        options::ComparisonOptions,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    qlen = length(query_variant.positions)
    slen = length(search_variant.positions)
    # Shift convention:
    # - positive shift moves search right relative to query,
    # - overlap indices computed in query coordinate space.
    overlap_start = max(1, 1 + shift)
    overlap_end = min(qlen, slen + shift)
    overlap_length = overlap_end - overlap_start + 1
    overlap_length < 1 && return nothing
    search_overlap_start = overlap_start - shift
    search_overlap_end = overlap_end - shift

    query_clips_fixed_prefix = _query_fixed_required(options.matchfix) &&
                               _clips_fixed_prefix(query_variant.positions, overlap_start)
    query_clips_fixed_suffix = _query_fixed_required(options.matchfix) &&
                               _clips_fixed_suffix(query_variant.positions, overlap_end)
    search_clips_fixed_prefix = _search_fixed_required(options.matchfix) &&
                                _clips_fixed_prefix(
        search_variant.positions,
        search_overlap_start
    )
    search_clips_fixed_suffix = _search_fixed_required(options.matchfix) &&
                                _clips_fixed_suffix(
        search_variant.positions,
        search_overlap_end
    )
    if (query_clips_fixed_prefix || query_clips_fixed_suffix ||
        search_clips_fixed_prefix || search_clips_fixed_suffix) &&
       !_is_exact_contained_alignment(
        query_variant.positions,
        search_variant.positions,
        overlap_start,
        search_overlap_start,
        overlap_length,
        query_clips_fixed_prefix,
        query_clips_fixed_suffix,
        search_clips_fixed_prefix,
        search_clips_fixed_suffix,
        spec.mask
    )
        return nothing
    end

    matched_pattern = IOBuffer()
    matched_positions = 0
    exact_fixed_matches = 0
    mismatches = 0
    match_ic = 0.0
    core_ic_denominator = 0.0
    has_variant = false
    has_degenerate = false
    has_complex = false

    for qidx in overlap_start:overlap_end
        sidx = qidx - shift
        qpos = query_variant.positions[qidx]
        spos = search_variant.positions[sidx]
        cmp = _compare_positions(qpos, spos, options, spec, residue_frequencies)
        if cmp.hard_mismatch
            return nothing
        end
        core_ic_denominator += cmp.core_ic_denominator

        if cmp.mismatch
            mismatches += 1
            mismatches > options.mismatches && return nothing
        else
            match_ic += cmp.ic
            matched_positions += cmp.contributes_position ? 1 : 0
            exact_fixed_matches += cmp.exact_fixed ? 1 : 0
            has_variant |= cmp.relation == _REL_VARIANT
            has_degenerate |= cmp.relation == _REL_DEGENERATE
            has_complex |= cmp.relation == _REL_COMPLEX
        end

        if _query_fixed_required(options.matchfix)
            if _is_fixed(qpos) && !cmp.exact_fixed
                # Required fixed query position must align to same fixed residue.
                return nothing
            end
        end
        if _search_fixed_required(options.matchfix)
            if _is_fixed(spos) && !cmp.exact_fixed
                # Required fixed search position must align to same fixed residue.
                return nothing
            end
        end

        print(matched_pattern,
            _match_symbol(
                qpos, spos, cmp.intersection, cmp.relation, cmp.mismatch, spec))
    end

    if matched_positions < options.min_shared_positions
        return nothing
    end

    # Normalize by the lower-information motif to keep score symmetric.
    denom = min(query_variant.information, search_variant.information)
    normalized_ic = denom > 0 ? (match_ic / denom) : 0.0
    if normalized_ic < options.normalized_ic_cutoff
        return nothing
    end

    # Oracle `CoreIC` is normalized by the more informative side of each aligned
    # position, not by a raw core-position count.
    core_ic = core_ic_denominator > 0 ? (match_ic / core_ic_denominator) : 0.0
    score = normalized_ic * matched_positions
    rel_type = _relationship_type_from_flags(has_variant, has_degenerate, has_complex)
    rel_length = _relationship_length(qlen, slen, overlap_length)
    rev_type = _reverse_type(rel_type)
    rev_length = _reverse_length(rel_length)

    return _Candidate(
        query_variant,
        search_variant,
        rel_type,
        rel_length,
        rev_type,
        rev_length,
        String(take!(matched_pattern)),
        matched_positions,
        exact_fixed_matches,
        match_ic,
        normalized_ic,
        core_ic,
        score
    )
end

"""
    _is_better(candidate::_Candidate, best::Union{Nothing, _Candidate})::Bool

Apply deterministic candidate ordering: 1) higher `match_ic`, 2) more matched positions,
3) higher `score`, 4) more exact fixed matches. Remaining ties fall back to
candidate encounter order from the shift scan inferred by black-box oracle tie cases.
"""
function _is_better(candidate::_Candidate, best::Union{Nothing, _Candidate})
    best === nothing && return true
    if candidate.match_ic != best.match_ic
        return candidate.match_ic > best.match_ic
    end
    if candidate.matched_positions != best.matched_positions
        return candidate.matched_positions > best.matched_positions
    end
    if candidate.score != best.score
        return candidate.score > best.score
    end
    if candidate.exact_fixed_matches != best.exact_fixed_matches
        return candidate.exact_fixed_matches > best.exact_fixed_matches
    end
    return false
end

"""
    _matches_exact_subsequence(shorter, longer, start)::Bool

Return `true` when `shorter` matches `longer[start:start+length(shorter)-1]`
position-for-position with identical anchors/residue classes.
"""
function _matches_exact_subsequence(
        shorter::Vector{_Position},
        longer::Vector{_Position},
        start::Int
)
    @inbounds for i in eachindex(shorter)
        apos = shorter[i]
        bpos = longer[start + i - 1]
        if apos.kind != bpos.kind || apos.mask != bpos.mask
            return false
        end
    end
    return true
end

"""
    _matches_matchfix_exact_subsequence(shorter, longer, start, wildcard_mask)::Bool

Return `true` when `shorter` is an exact contained subsequence of `longer` and
the shorter motif keeps an informative edge residue on every side where the
constrained longer motif clips fixed positions. Wildcards away from those
clipped boundaries are still allowed by the oracle.
"""
function _matches_matchfix_exact_subsequence(
        shorter::AbstractVector{_Position},
        longer::AbstractVector{_Position},
        start::Int,
        wildcard_mask::ResidueMask,
        require_leading_informative::Bool,
        require_trailing_informative::Bool
)
    @inbounds for i in eachindex(shorter)
        apos = shorter[i]
        bpos = longer[start + i - 1]
        if apos.kind != bpos.kind || apos.mask != bpos.mask
            return false
        end
    end
    require_leading_informative &&
        _is_wildcard(shorter[firstindex(shorter)], wildcard_mask) && return false
    require_trailing_informative &&
        _is_wildcard(shorter[lastindex(shorter)], wildcard_mask) && return false
    return true
end

@inline function _first_informative_index(
        positions::AbstractVector{_Position},
        wildcard_mask::ResidueMask
)
    @inbounds for idx in eachindex(positions)
        !_is_wildcard(positions[idx], wildcard_mask) && return idx
    end
    return nothing
end

@inline function _last_informative_index(
        positions::AbstractVector{_Position},
        wildcard_mask::ResidueMask
)
    @inbounds for idx in reverse(eachindex(positions))
        !_is_wildcard(positions[idx], wildcard_mask) && return idx
    end
    return nothing
end

function _clipped_region_can_rebind(
        positions::AbstractVector{_Position},
        start::Int,
        stop::Int,
        target::_Position
)
    start > stop && return false
    @inbounds for idx in start:stop
        pos = positions[idx]
        if pos.kind != target.kind
            continue
        end
        if pos.kind != _RESIDUE ||
           overlaps(ResidueClass(pos.mask), ResidueClass(target.mask))
            return true
        end
    end
    return false
end

"""
    _matches_matchfix_precise_subsequence(shorter, longer, start, longer_constrained, wildcard_mask)::Bool

Return `true` when an exact contained alignment should still be classified as an
exact subsequence/parent under `matchfix`. When a wildcard sits on the clipped
edge of the shorter motif, the oracle only keeps the exact classification if
either the longer side is unconstrained and the clipped residues cannot rebind
the nearest informative edge position, or the wildcard is on the opposite edge.
"""
function _matches_matchfix_precise_subsequence(
        shorter::AbstractVector{_Position},
        longer::AbstractVector{_Position},
        start::Int,
        longer_constrained::Bool,
        wildcard_mask::ResidueMask
)
    _matches_exact_subsequence(shorter, longer, start) || return false

    prefix_clipped = start > firstindex(longer)
    suffix_clipped = (start + length(shorter) - 1) < lastindex(longer)

    if prefix_clipped && _is_wildcard(shorter[firstindex(shorter)], wildcard_mask)
        longer_constrained && return false
        informative_idx = _first_informative_index(shorter, wildcard_mask)
        informative_idx === nothing && return false
        _clipped_region_can_rebind(
            longer,
            firstindex(longer),
            start - 1,
            shorter[informative_idx]
        ) && return false
    end

    if suffix_clipped && _is_wildcard(shorter[lastindex(shorter)], wildcard_mask)
        longer_constrained && return false
        informative_idx = _last_informative_index(shorter, wildcard_mask)
        informative_idx === nothing && return false
        _clipped_region_can_rebind(
            longer,
            start + length(shorter),
            lastindex(longer),
            shorter[informative_idx]
        ) && return false
    end

    return true
end

"""
    _find_precise_match(query_variants, search_variants, options, spec)

Search only exact same / exact-subsequence relationships. Returns
`(found_precise, best_candidate)`.
"""
function _find_precise_match(
        query_variants::Vector{_MotifVariant},
        search_variants::Vector{_MotifVariant},
        options::ComparisonOptions,
        spec::_AlphabetSpec
)
    return _find_precise_match(
        query_variants,
        search_variants,
        options,
        spec,
        _residue_frequency_vector(options, spec)
    )
end

function _find_precise_match(
        query_variants::Vector{_MotifVariant},
        search_variants::Vector{_MotifVariant},
        options::ComparisonOptions,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    found_precise = false
    best = nothing

    for qvariant in query_variants
        qlen = length(qvariant.positions)
        for svariant in search_variants
            slen = length(svariant.positions)

            if qlen <= slen
                for start in 1:(slen - qlen + 1)
                    if !_matches_matchfix_precise_subsequence(
                        qvariant.positions,
                        svariant.positions,
                        start,
                        _search_fixed_required(options.matchfix),
                        spec.mask
                    )
                        continue
                    end
                    found_precise = true
                    candidate = _evaluate_alignment(
                        qvariant, svariant, 1 - start, options, spec, residue_frequencies)
                    if candidate !== nothing && _is_better(candidate, best)
                        best = candidate
                    end
                end
            end

            if slen < qlen
                for start in 1:(qlen - slen + 1)
                    if !_matches_matchfix_precise_subsequence(
                        svariant.positions,
                        qvariant.positions,
                        start,
                        _query_fixed_required(options.matchfix),
                        spec.mask
                    )
                        continue
                    end
                    found_precise = true
                    candidate = _evaluate_alignment(
                        qvariant, svariant, start - 1, options, spec, residue_frequencies)
                    if candidate !== nothing && _is_better(candidate, best)
                        best = candidate
                    end
                end
            end
        end
    end

    return (found_precise, best)
end

"""
    _compare_parsed(parsed_query, parsed_search, options)::ComparisonResult

Compare two already-parsed motifs.
"""
function _compare_parsed(parsed_query::_ParsedMotif, parsed_search::_ParsedMotif, options::ComparisonOptions)
    # Repeat ranges are expanded first; alignment runs over concrete variants.
    spec = _alphabet_spec(options.alphabet)
    residue_frequencies = _residue_frequency_vector(options, spec)
    query_variants = _expand_variants(parsed_query, options, spec, residue_frequencies)
    search_variants = _expand_variants(parsed_search, options, spec, residue_frequencies)
    _,
    best = _find_precise_match(
        query_variants,
        search_variants,
        options,
        spec,
        residue_frequencies
    )

    # The paper's "precise match first" rule is applied after regex
    # preprocessing/expansion, but the best hit is still selected across the
    # full expanded motif pair. Exact-subsequence hits from one branch must not
    # suppress stronger overlaps from other branches.
    for qvariant in query_variants
        for svariant in search_variants
            qlen = length(qvariant.positions)
            slen = length(svariant.positions)
            for shift in (qlen - 1):-1:(-(slen - 1))
                candidate = _evaluate_alignment(
                    qvariant, svariant, shift, options, spec, residue_frequencies)
                candidate === nothing && continue
                if _is_better(candidate, best)
                    best = candidate
                end
            end
        end
    end

    if best === nothing
        return ComparisonResult(
            query = parsed_query.original,
            search = parsed_search.original,
            normalized_query = parsed_query.normalized,
            normalized_search = parsed_search.normalized
        )
    end

    # Materialize public result fields from the best candidate.
    return ComparisonResult(
        query = parsed_query.original,
        search = parsed_search.original,
        normalized_query = parsed_query.normalized,
        normalized_search = parsed_search.normalized,
        matched = true,
        query_relationship = _relationship_word(
            best.query_relationship_type, best.query_relationship_length),
        search_relationship = _relationship_word(
            best.search_relationship_type, best.search_relationship_length),
        matched_pattern = best.matched_pattern,
        matched_positions = best.matched_positions,
        match_ic = best.match_ic,
        normalized_ic = best.normalized_ic,
        core_ic = best.core_ic,
        score = best.score,
        query_information = best.query_variant.information,
        search_information = best.search_variant.information
    )
end
