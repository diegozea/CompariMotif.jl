"""
    _anchor_symbol(pos::_Position)::String

Render the canonical motif symbol for an anchor position (`^` or `\$`).
"""
@inline function _anchor_symbol(pos::_Position)
    pos.kind == _NTERMINUS && return "^"
    pos.kind == _CTERMINUS && return "\$"
    error("Expected anchor position")
end

"""
    _anchor_mismatch_symbol(anchor_pos::_Position, residue_pos::_Position, spec::_AlphabetSpec)::String

Render the overlap symbol used when an anchor aligns against a residue-bearing
position.
"""
function _anchor_mismatch_symbol(
        anchor_pos::_Position,
        residue_pos::_Position,
        spec::_AlphabetSpec
)
    anchor = _anchor_symbol(anchor_pos)
    if _is_wildcard(residue_pos, spec.mask)
        return anchor
    end
    residues = String(_mask_to_chars(residue_pos.mask, spec; as_lowercase = true))
    return "[$anchor$residues]"
end

"""
    _position_comparison(; kwargs...)

Build the named-tuple returned by position-comparison helpers.

The fields encode whether the aligned position is a hard rejection, whether it
consumes mismatch budget, the per-position relationship label, and the scoring
terms accumulated by `_evaluate_alignment`.
"""
@inline function _position_comparison(;
        hard_mismatch::Bool,
        mismatch::Bool,
        relation::_RelationshipType,
        intersection::ResidueMask = ResidueMask(0),
        ic::Float64 = 0.0,
        core_ic_denominator::Float64 = 0.0,
        contributes_position::Bool = false,
        exact_fixed::Bool = false
)
    return (
        hard_mismatch = hard_mismatch,
        mismatch = mismatch,
        relation = relation,
        intersection = intersection,
        ic = ic,
        core_ic_denominator = core_ic_denominator,
        contributes_position = contributes_position,
        exact_fixed = exact_fixed
    )
end

"""
    _compare_anchor_positions(qpos::_Position, spos::_Position)

Compare two anchor positions. Matching anchors produce an exact informative
position; opposing anchors are treated as a hard mismatch.
"""
function _compare_anchor_positions(qpos::_Position, spos::_Position)
    if qpos.kind == spos.kind
        return _position_comparison(;
            hard_mismatch = false,
            mismatch = false,
            relation = _REL_EXACT,
            ic = 1.0,
            core_ic_denominator = 1.0,
            contributes_position = true
        )
    end
    return _position_comparison(;
        hard_mismatch = true,
        mismatch = true,
        relation = _REL_COMPLEX
    )
end

"""
    _compare_anchor_residue_positions(qpos, spos, options, spec, residue_frequencies)

Compare an anchor against a residue-bearing position.

Under mismatch-tolerant mode this usually consumes one mismatch while still
contributing the larger anchor/residue information term to the CoreIC
denominator. When the anchored side is constrained by `matchfix`, the alignment
is rejected outright.
"""
function _compare_anchor_residue_positions(
        qpos::_Position,
        spos::_Position,
        options::ComparisonOptions,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    constrained_anchor_mismatch = (qpos.kind !== _RESIDUE &&
                                   _query_fixed_required(options.matchfix)) ||
                                  (spos.kind !== _RESIDUE &&
                                   _search_fixed_required(options.matchfix))
    constrained_anchor_mismatch && return _position_comparison(;
        hard_mismatch = true,
        mismatch = true,
        relation = _REL_COMPLEX
    )

    q_ic = _position_ic(qpos, spec, residue_frequencies)
    s_ic = _position_ic(spos, spec, residue_frequencies)
    return _position_comparison(;
        hard_mismatch = false,
        mismatch = true,
        relation = _REL_COMPLEX,
        core_ic_denominator = max(q_ic, s_ic)
    )
end

"""
    _position_relation(qclass, sclass, allow_ambiguous_overlap)

Classify the residue-set relationship for one aligned residue position.

Returns one of the internal relationship tags, or `nothing` when the classes
only partially overlap and ambiguous overlaps are disallowed.
"""
function _position_relation(
        qclass::ResidueClass,
        sclass::ResidueClass,
        allow_ambiguous_overlap::Bool
)
    if qclass.mask == sclass.mask
        return _REL_EXACT
    elseif is_subset(qclass, sclass)
        return _REL_VARIANT
    elseif is_subset(sclass, qclass)
        return _REL_DEGENERATE
    elseif allow_ambiguous_overlap
        return _REL_COMPLEX
    end
    return nothing
end

"""
    _position_match_ic(relation, qclass, sclass, q_ic, s_ic, spec, residue_frequencies)::Float64

Return the information content credited to one matched residue position.

Exact/variant/degenerate matches use the less informative side; complex partial
overlaps use the union-class information content to match the oracle's scoring.
"""
function _position_match_ic(
        relation::_RelationshipType,
        qclass::ResidueClass,
        sclass::ResidueClass,
        q_ic::Float64,
        s_ic::Float64,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    if relation == _REL_COMPLEX
        union_pos = _Position(_RESIDUE, unionclass(qclass, sclass).mask)
        return _position_ic(union_pos, spec, residue_frequencies)
    end
    return min(q_ic, s_ic)
end

"""
    _compare_residue_positions(qpos, spos, options, spec, residue_frequencies)

Compare two residue-bearing positions and return the diagnostic tuple consumed
by `_evaluate_alignment`.
"""
function _compare_residue_positions(
        qpos::_Position,
        spos::_Position,
        options::ComparisonOptions,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    q_ic = _position_ic(qpos, spec, residue_frequencies)
    s_ic = _position_ic(spos, spec, residue_frequencies)

    qclass = ResidueClass(qpos.mask)
    sclass = ResidueClass(spos.mask)
    intersection = qclass.mask & sclass.mask
    core_ic_denominator = max(q_ic, s_ic)

    overlaps(qclass, sclass) || return _position_comparison(;
        hard_mismatch = false,
        mismatch = true,
        relation = _REL_COMPLEX,
        intersection = intersection,
        core_ic_denominator = core_ic_denominator
    )

    relation = _position_relation(qclass, sclass, options.allow_ambiguous_overlap)
    relation === nothing && return _position_comparison(;
        hard_mismatch = true,
        mismatch = true,
        relation = _REL_COMPLEX,
        intersection = intersection,
        core_ic_denominator = core_ic_denominator
    )

    ic = _position_match_ic(relation, qclass, sclass, q_ic, s_ic, spec, residue_frequencies)
    contributes_position = !_is_wildcard(qpos, spec.mask) && !_is_wildcard(spos, spec.mask)
    exact_fixed = relation == _REL_EXACT && _is_fixed(qpos) && _is_fixed(spos)
    return _position_comparison(;
        hard_mismatch = false,
        mismatch = false,
        relation = relation,
        intersection = intersection,
        ic = ic,
        core_ic_denominator = core_ic_denominator,
        contributes_position = contributes_position,
        exact_fixed = exact_fixed
    )
end

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
    if qpos.kind !== _RESIDUE && spos.kind !== _RESIDUE
        if qpos.kind == spos.kind
            return _anchor_symbol(qpos)
        end
        return "."
    elseif qpos.kind !== _RESIDUE
        return _anchor_mismatch_symbol(qpos, spos, spec)
    elseif spos.kind !== _RESIDUE
        return _anchor_mismatch_symbol(spos, qpos, spec)
    end

    qwild = _is_wildcard(qpos, spec.mask)
    swild = _is_wildcard(spos, spec.mask)
    if qwild && swild
        # Preserve compact wildcard representation in overlap output.
        return "."
    end

    qclass = ResidueClass(qpos.mask)
    sclass = ResidueClass(spos.mask)
    if mismatch
        # Mismatch still contributes a descriptive symbol in the matched pattern.
        union_mask = unionclass(qclass, sclass).mask
        return _mask_to_symbol(union_mask, spec; as_lowercase = true)
    end
    if relation == _REL_EXACT
        return _mask_to_symbol(intersection, spec)
    end
    if qwild
        return _mask_to_symbol(spos.mask, spec; as_lowercase = true)
    elseif swild
        return _mask_to_symbol(qpos.mask, spec; as_lowercase = true)
    end
    union_mask = unionclass(qclass, sclass).mask
    return _mask_to_symbol(union_mask, spec; as_lowercase = true)
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
    if qpos.kind !== _RESIDUE && spos.kind !== _RESIDUE
        return _compare_anchor_positions(qpos, spos)
    elseif qpos.kind !== _RESIDUE || spos.kind !== _RESIDUE
        return _compare_anchor_residue_positions(
            qpos, spos, options, spec, residue_frequencies)
    end
    return _compare_residue_positions(qpos, spos, options, spec, residue_frequencies)
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

# In `matchfix` modes, the constrained side cannot hide informative positions
# outside the overlap. "Fixed" means a single defined residue and "clipped"
# means that the chosen alignment window leaves such a position on the prefix or
# suffix outside the aligned span.
"""
    _clips_fixed_prefix(positions, overlap_start)::Bool

Return `true` when the constrained side clips one or more fixed residues before
the aligned span.

```jldoctest
julia> using CompariMotif

julia> options = CompariMotif.ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> spec = CompariMotif._alphabet_spec(options.alphabet);

julia> variant = only(CompariMotif._expand_variants(
           CompariMotif._parse_motif(raw"^ACD\$", options),
           options,
           spec
       ));

julia> CompariMotif._clips_fixed_prefix(variant.positions, 3)
true
```
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
    _clips_anchor_prefix(positions, overlap_start)::Bool

Return `true` when the constrained side clips one or more anchors before the
aligned span.
"""
function _clips_anchor_prefix(
        positions::AbstractVector{_Position},
        overlap_start::Int
)
    @inbounds for idx in firstindex(positions):(overlap_start - 1)
        if positions[idx].kind !== _RESIDUE
            return true
        end
    end
    return false
end

"""
    _clips_anchor_suffix(positions, overlap_end)::Bool

Return `true` when the constrained side clips one or more anchors after the
aligned span.
"""
function _clips_anchor_suffix(
        positions::AbstractVector{_Position},
        overlap_end::Int
)
    @inbounds for idx in (overlap_end + 1):lastindex(positions)
        if positions[idx].kind !== _RESIDUE
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
    _alignment_window(query_positions, search_positions, shift)

Compute the overlap geometry for one query/search shift in query coordinates.

Returns `nothing` when the shift produces no overlap; otherwise returns a named
tuple containing both query- and search-space overlap bounds.
"""
function _alignment_window(
        query_positions::AbstractVector{_Position},
        search_positions::AbstractVector{_Position},
        shift::Int
)
    qlen = length(query_positions)
    slen = length(search_positions)
    overlap_start = max(1, 1 + shift)
    overlap_end = min(qlen, slen + shift)
    overlap_length = overlap_end - overlap_start + 1
    overlap_length < 1 && return nothing
    return (
        qlen = qlen,
        slen = slen,
        overlap_start = overlap_start,
        overlap_end = overlap_end,
        overlap_length = overlap_length,
        search_overlap_start = overlap_start - shift,
        search_overlap_end = overlap_end - shift
    )
end

"""
    _alignment_clip_state(positions, overlap_start, overlap_end, constrained)

Summarize whether a constrained motif clips fixed residues or anchors on either
side of the chosen overlap window.

When `constrained == false`, all clip flags are `false`.
"""
function _alignment_clip_state(
        positions::AbstractVector{_Position},
        overlap_start::Int,
        overlap_end::Int,
        constrained::Bool
)
    constrained || return (
        fixed_prefix = false,
        fixed_suffix = false,
        anchor_prefix = false,
        anchor_suffix = false
    )
    return (
        fixed_prefix = _clips_fixed_prefix(positions, overlap_start),
        fixed_suffix = _clips_fixed_suffix(positions, overlap_end),
        anchor_prefix = _clips_anchor_prefix(positions, overlap_start),
        anchor_suffix = _clips_anchor_suffix(positions, overlap_end)
    )
end

"""
    _clips_fixed(state)::Bool

Return `true` when either side of a clip-state summary hides fixed residues
outside the aligned span.
"""
@inline _clips_fixed(state) = state.fixed_prefix || state.fixed_suffix

"""
    _clips_anchor(state)::Bool

Return `true` when either side of a clip-state summary hides anchors outside
the aligned span.
"""
@inline _clips_anchor(state) = state.anchor_prefix || state.anchor_suffix

"""
    _rejects_clipped_fixed_alignment(query_positions, search_positions, overlap, query_clip_state, search_clip_state, wildcard_mask)::Bool

Return `true` when `matchfix` clipped-fixed rules reject the current overlap
before per-position scoring begins.
"""
function _rejects_clipped_fixed_alignment(
        query_positions::AbstractVector{_Position},
        search_positions::AbstractVector{_Position},
        overlap,
        query_clip_state,
        search_clip_state,
        wildcard_mask::ResidueMask
)
    (_clips_fixed(query_clip_state) || _clips_fixed(search_clip_state)) || return false
    return !_is_exact_contained_alignment(
        query_positions,
        search_positions,
        overlap.overlap_start,
        overlap.search_overlap_start,
        overlap.overlap_length,
        query_clip_state.fixed_prefix,
        query_clip_state.fixed_suffix,
        search_clip_state.fixed_prefix,
        search_clip_state.fixed_suffix,
        wildcard_mask
    )
end

"""
    _AlignmentAccumulator

Mutable aggregation state for one candidate overlap while `_evaluate_alignment`
scans positions from left to right.

Fields:
- `matched_pattern`: overlap text being emitted for the winning span.
- `matched_positions`: informative non-wildcard matched-position count.
- `mismatches`: mismatches consumed so far.
- `match_ic`: accumulated matched information content.
- `core_ic_denominator`: accumulated CoreIC denominator.
- `has_variant`, `has_degenerate`, `has_complex`: per-position relationship
  evidence used to derive the final relationship word.
"""
Base.@kwdef mutable struct _AlignmentAccumulator
    matched_pattern::IOBuffer = IOBuffer()
    matched_positions::Int = 0
    mismatches::Int = 0
    match_ic::Float64 = 0.0
    core_ic_denominator::Float64 = 0.0
    has_variant::Bool = false
    has_degenerate::Bool = false
    has_complex::Bool = false
end

"""
    _rejects_clipped_anchor_mismatch(qpos, spos, query_clip_state, search_clip_state)::Bool

Return `true` when an anchor-vs-residue alignment should be rejected because a
constrained side already clips anchors elsewhere in the overlap.
"""
@inline function _rejects_clipped_anchor_mismatch(
        qpos::_Position,
        spos::_Position,
        query_clip_state,
        search_clip_state
)
    ((qpos.kind !== _RESIDUE) ⊻ (spos.kind !== _RESIDUE)) || return false
    return _clips_anchor(query_clip_state) || _clips_anchor(search_clip_state)
end

"""
    _violates_matchfix_exact_fixed(qpos, spos, cmp, query_constrained, search_constrained)::Bool

Return `true` when a fixed residue on a constrained side fails to align to the
same fixed residue exactly.
"""
@inline function _violates_matchfix_exact_fixed(
        qpos::_Position,
        spos::_Position,
        cmp,
        query_constrained::Bool,
        search_constrained::Bool
)
    return (query_constrained && _is_fixed(qpos) && !cmp.exact_fixed) ||
           (search_constrained && _is_fixed(spos) && !cmp.exact_fixed)
end

"""
    _record_alignment_position!(acc, qpos, spos, cmp, options, spec)::Bool

Update the running overlap accumulator with one aligned position.

Returns `false` only when recording the position would exceed the allowed
mismatch budget.
"""
function _record_alignment_position!(
        acc::_AlignmentAccumulator,
        qpos::_Position,
        spos::_Position,
        cmp,
        options::ComparisonOptions,
        spec::_AlphabetSpec
)
    acc.core_ic_denominator += cmp.core_ic_denominator
    if cmp.mismatch
        acc.mismatches += 1
        acc.mismatches > options.mismatches && return false
    else
        acc.match_ic += cmp.ic
        acc.matched_positions += cmp.contributes_position ? 1 : 0
        acc.has_variant |= cmp.relation == _REL_VARIANT
        acc.has_degenerate |= cmp.relation == _REL_DEGENERATE
        acc.has_complex |= cmp.relation == _REL_COMPLEX
    end

    print(acc.matched_pattern,
        _match_symbol(qpos, spos, cmp.intersection, cmp.relation, cmp.mismatch, spec))
    return true
end

"""
    _build_alignment_candidate(query_variant, search_variant, overlap_length, acc, options)

Finalize one accumulated overlap into a scored `_Candidate`.

Returns `nothing` when the overlap fails the global matched-position or
normalized-IC thresholds.
"""
function _build_alignment_candidate(
        query_variant::_MotifVariant,
        search_variant::_MotifVariant,
        overlap_length::Int,
        acc::_AlignmentAccumulator,
        options::ComparisonOptions
)
    acc.matched_positions < options.min_shared_positions && return nothing

    denom = min(query_variant.information, search_variant.information)
    normalized_ic = denom > 0 ? (acc.match_ic / denom) : 0.0
    normalized_ic < options.normalized_ic_cutoff && return nothing

    core_ic = acc.core_ic_denominator > 0 ? (acc.match_ic / acc.core_ic_denominator) : 0.0
    score = normalized_ic * acc.matched_positions
    qlen = length(query_variant.positions)
    slen = length(search_variant.positions)
    rel_type = _relationship_type_from_flags(
        acc.has_variant, acc.has_degenerate, acc.has_complex)
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
        String(take!(acc.matched_pattern)),
        acc.matched_positions,
        acc.match_ic,
        normalized_ic,
        core_ic,
        score
    )
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
    overlap = _alignment_window(query_variant.positions, search_variant.positions, shift)
    overlap === nothing && return nothing

    query_constrained = _query_fixed_required(options.matchfix)
    search_constrained = _search_fixed_required(options.matchfix)
    query_clip_state = _alignment_clip_state(
        query_variant.positions,
        overlap.overlap_start,
        overlap.overlap_end,
        query_constrained
    )
    search_clip_state = _alignment_clip_state(
        search_variant.positions,
        overlap.search_overlap_start,
        overlap.search_overlap_end,
        search_constrained
    )
    _rejects_clipped_fixed_alignment(
        query_variant.positions,
        search_variant.positions,
        overlap,
        query_clip_state,
        search_clip_state,
        spec.mask
    ) && return nothing

    acc = _AlignmentAccumulator()
    for qidx in overlap.overlap_start:overlap.overlap_end
        sidx = qidx - shift
        qpos = query_variant.positions[qidx]
        spos = search_variant.positions[sidx]
        cmp = _compare_positions(qpos, spos, options, spec, residue_frequencies)
        cmp.hard_mismatch && return nothing
        _rejects_clipped_anchor_mismatch(qpos, spos, query_clip_state, search_clip_state) &&
            return nothing
        _violates_matchfix_exact_fixed(
            qpos, spos, cmp, query_constrained, search_constrained) &&
            return nothing
        _record_alignment_position!(acc, qpos, spos, cmp, options, spec) || return nothing
    end

    return _build_alignment_candidate(
        query_variant,
        search_variant,
        overlap.overlap_length,
        acc,
        options
    )
end

"""
    _is_better(candidate::_Candidate, best::Union{Nothing, _Candidate})::Bool

Apply deterministic candidate ordering: 1) higher `match_ic`, 2) more matched positions,
3) higher `score`. Remaining ties fall back to candidate encounter order from
the shift scan inferred by black-box oracle tie cases.
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

"""
    _first_informative_index(positions, wildcard_mask)

Return the first index whose position is not a wildcard, or `nothing` when the
entire sequence is wildcard-only.
"""
@inline function _first_informative_index(
        positions::AbstractVector{_Position},
        wildcard_mask::ResidueMask
)
    @inbounds for idx in eachindex(positions)
        !_is_wildcard(positions[idx], wildcard_mask) && return idx
    end
    return nothing
end

"""
    _last_informative_index(positions, wildcard_mask)

Return the last index whose position is not a wildcard, or `nothing` when the
entire sequence is wildcard-only.
"""
@inline function _last_informative_index(
        positions::AbstractVector{_Position},
        wildcard_mask::ResidueMask
)
    @inbounds for idx in reverse(eachindex(positions))
        !_is_wildcard(positions[idx], wildcard_mask) && return idx
    end
    return nothing
end

"""
    _clipped_region_can_rebind(positions, start, stop, target) :: Bool

Return `true` when a clipped region contains an anchor of the same kind or a
residue class overlapping `target`, meaning the clipped segment could compete
for that edge position under the oracle's precise-match rules.
"""
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
