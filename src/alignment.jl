"""
    _match_symbol(qpos, spos, intersection, relation, mismatch, options) -> String

Render one output symbol for the overlap pattern.
"""
function _match_symbol(
        qpos::_Position,
        spos::_Position,
        intersection::ResidueMask,
        relation::_RelationshipType,
        mismatch::Bool,
        options::ComparisonOptions
)
    if qpos.kind == _NTERMINUS
        return "^"
    elseif qpos.kind == _CTERMINUS
        return "\$"
    end

    qwild = _is_wildcard(qpos, options)
    swild = _is_wildcard(spos, options)
    if qwild && swild
        # Preserve compact wildcard representation in overlap output.
        return "x"
    end

    qclass = ResidueClass(qpos.mask)
    sclass = ResidueClass(spos.mask)
    if mismatch
        # Mismatch still contributes a descriptive symbol in the matched pattern.
        union_mask = unionclass(qclass, sclass).mask
        return _mask_to_symbol(union_mask, options; as_lowercase = true, wildcard_symbol = "x")
    end
    if relation == _REL_EXACT
        if qwild || swild
            if qwild
                return _mask_to_symbol(spos.mask, options; as_lowercase = true, wildcard_symbol = "x")
            end
            return _mask_to_symbol(qpos.mask, options; as_lowercase = true, wildcard_symbol = "x")
        end
        return _mask_to_symbol(intersection, options; as_lowercase = false, wildcard_symbol = "x")
    end
    if qwild
        return _mask_to_symbol(spos.mask, options; as_lowercase = true, wildcard_symbol = "x")
    elseif swild
        return _mask_to_symbol(qpos.mask, options; as_lowercase = true, wildcard_symbol = "x")
    end
    union_mask = unionclass(qclass, sclass).mask
    return _mask_to_symbol(union_mask, options; as_lowercase = true, wildcard_symbol = "x")
end

"""
    _compare_positions(qpos, spos, options)

Compare one query/search position pair and return matching diagnostics.
"""
function _compare_positions(qpos::_Position, spos::_Position, options::ComparisonOptions)
    # Anchors are valid only against the same anchor type.
    if qpos.kind !== _RESIDUE || spos.kind !== _RESIDUE
        if qpos.kind == spos.kind
            return (hard_mismatch = false, mismatch = false,
                relation = _REL_EXACT, intersection = ResidueMask(0), ic = 1.0,
                contributes_position = true, exact_fixed = false)
        end
        return (hard_mismatch = true, mismatch = true,
            relation = _REL_COMPLEX, intersection = ResidueMask(0), ic = 0.0,
            contributes_position = false, exact_fixed = false)
    end

    qclass = ResidueClass(qpos.mask)
    sclass = ResidueClass(spos.mask)
    # Residue match operates on mask set intersection.
    intersection = qclass.mask & sclass.mask
    if !overlaps(qclass, sclass)
        # No shared residue -> mismatch position (allowed only if mismatch budget allows).
        return (hard_mismatch = false, mismatch = true,
            relation = _REL_COMPLEX, intersection = intersection, ic = 0.0,
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
            contributes_position = false, exact_fixed = false)
    end

    # Position IC contribution follows the less-specific side.
    ic = min(_position_ic(qpos, options), _position_ic(spos, options))
    # Matched position count excludes wildcard-vs-wildcard.
    contributes = !_is_wildcard(qpos, options) && !_is_wildcard(spos, options)
    # Used for matchfix tie-break logic.
    exact_fixed = relation == _REL_EXACT && _is_fixed(qpos) && _is_fixed(spos)
    return (
        hard_mismatch = false,
        mismatch = false,
        relation = relation,
        intersection = intersection,
        ic = ic,
        contributes_position = contributes,
        exact_fixed = exact_fixed
    )
end

"""
    _query_fixed_required(mode::MatchFixMode) -> Bool

Return `true` when query fixed residues must match exactly.
"""
function _query_fixed_required(mode::MatchFixMode)
    # Query fixed constraints apply in QueryFixed and BothFixed modes.
    mode == MatchFixQueryFixed || mode == MatchFixBothFixed
end

"""
    _search_fixed_required(mode::MatchFixMode) -> Bool

Return `true` when search fixed residues must match exactly.
"""
function _search_fixed_required(mode::MatchFixMode)
    # Search fixed constraints apply in SearchFixed and BothFixed modes.
    mode == MatchFixSearchFixed || mode == MatchFixBothFixed
end

"""
    _evaluate_alignment(query_variant, search_variant, shift, options)

Evaluate one concrete shift between two expanded motif variants.
Returns `_Candidate` when all thresholds pass, otherwise `nothing`.
"""
function _evaluate_alignment(
        query_variant::_MotifVariant,
        search_variant::_MotifVariant,
        shift::Int,
        options::ComparisonOptions
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

    matched_pattern = IOBuffer()
    matched_positions = 0
    exact_fixed_matches = 0
    mismatches = 0
    match_ic = 0.0
    core_length = 0
    has_variant = false
    has_degenerate = false
    has_complex = false

    for qidx in overlap_start:overlap_end
        sidx = qidx - shift
        qpos = query_variant.positions[qidx]
        spos = search_variant.positions[sidx]
        cmp = _compare_positions(qpos, spos, options)
        if cmp.hard_mismatch
            return nothing
        end
        if !(_is_wildcard(qpos, options) && _is_wildcard(spos, options))
            # Core length ignores dual-wildcard positions.
            core_length += 1
        end

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
                qpos, spos, cmp.intersection, cmp.relation, cmp.mismatch, options))
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

    core_ic = core_length > 0 ? (match_ic / core_length) : 0.0
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
    _is_better(candidate::_Candidate, best::Union{Nothing, _Candidate}) -> Bool

Apply deterministic candidate ordering:
1) higher `match_ic`, 2) more matched positions, 3) more exact fixed matches.
"""
function _is_better(candidate::_Candidate, best::Union{Nothing, _Candidate})
    best === nothing && return true
    if candidate.match_ic != best.match_ic
        return candidate.match_ic > best.match_ic
    end
    if candidate.matched_positions != best.matched_positions
        return candidate.matched_positions > best.matched_positions
    end
    if candidate.exact_fixed_matches != best.exact_fixed_matches
        return candidate.exact_fixed_matches > best.exact_fixed_matches
    end
    return false
end

"""
    _compare_parsed(parsed_query, parsed_search, options) -> ComparisonResult

Compare two already-parsed motifs.
"""
function _compare_parsed(parsed_query::_ParsedMotif, parsed_search::_ParsedMotif, options::ComparisonOptions)
    # Repeat ranges are expanded first; alignment runs over concrete variants.
    query_variants = _expand_variants(parsed_query, options)
    search_variants = _expand_variants(parsed_search, options)
    best = nothing

    # Evaluate all relative shifts for each variant pair.
    for qvariant in query_variants
        for svariant in search_variants
            qlen = length(qvariant.positions)
            slen = length(svariant.positions)
            for shift in (-(slen - 1)):(qlen - 1)
                candidate = _evaluate_alignment(qvariant, svariant, shift, options)
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
