_coerce_matchfix(mode::MatchFixMode) = mode
function _coerce_matchfix(mode::Symbol)
    # Accept a few user-friendly aliases while keeping one canonical enum internally.
    if mode === :none
        return MatchFixNone
    elseif mode === :query_fixed || mode === :query
        return MatchFixQueryFixed
    elseif mode === :search_fixed || mode === :search
        return MatchFixSearchFixed
    elseif mode === :both_fixed || mode === :both
        return MatchFixBothFixed
    end
    throw(ArgumentError("`matchfix` must be one of :none, :query_fixed, :search_fixed, :both_fixed."))
end
function _coerce_matchfix(mode::AbstractString)
    # Normalize spelling/case before reusing Symbol-based coercion.
    _coerce_matchfix(Symbol(replace(lowercase(strip(mode)), ' ' => '_')))
end

function ComparisonOptions(;
        alphabet::Symbol = :protein,
        min_shared_positions::Int = 2,
        normalized_ic_cutoff::Real = 0.5,
        matchfix::Union{MatchFixMode, Symbol, AbstractString} = MatchFixNone,
        mismatches::Int = 0,
        allow_ambiguous_overlap::Bool = true,
        max_variants::Int = 10_000
)
    # Validate and canonicalize simple scalar options first.
    matchfix_mode = _coerce_matchfix(matchfix)
    if min_shared_positions < 1
        throw(ArgumentError("`min_shared_positions` must be >= 1."))
    end
    if normalized_ic_cutoff < 0
        throw(ArgumentError("`normalized_ic_cutoff` must be >= 0."))
    end
    if mismatches < 0
        throw(ArgumentError("`mismatches` must be >= 0."))
    end
    if max_variants < 1
        throw(ArgumentError("`max_variants` must be >= 1."))
    end

    # Choose alphabet once and use it for both parser masks and IC formulas.
    alphabet_chars = if alphabet === :protein
        _PROTEIN_ALPHABET
    elseif alphabet === :dna
        _DNA_ALPHABET
    else
        throw(ArgumentError("`alphabet` must be :protein or :dna."))
    end

    # Build 1-based residue index -> bit index map.
    index = Dict{Char, Int}()
    for (i, aa) in enumerate(alphabet_chars)
        index[aa] = i
    end
    # Full alphabet mask: every residue bit set to one.
    mask = UInt32((1 << length(alphabet_chars)) - 1)

    return ComparisonOptions(
        alphabet_chars,
        index,
        mask,
        log(length(alphabet_chars)),
        min_shared_positions,
        float(normalized_ic_cutoff),
        matchfix_mode,
        mismatches,
        allow_ambiguous_overlap,
        max_variants
    )
end
