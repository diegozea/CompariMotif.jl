module CompariMotif

using AutoPrettyPrinting

export ComparisonOptions, ComparisonResult, MatchFixMode, MatchFixBothFixed,
       MatchFixNone, MatchFixQueryFixed, MatchFixSearchFixed, compare,
       normalize_motif, write_results_tsv

const _PROTEIN_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")
const _DNA_ALPHABET = collect("ACGT")

Base.@kwdef struct ComparisonResult
    query::String
    search::String
    normalized_query::String
    normalized_search::String
    matched::Bool = false
    query_relationship::String = "No Match"
    search_relationship::String = "No Match"
    query_relationship_code::String = ""
    search_relationship_code::String = ""
    matched_pattern::String = ""
    matched_positions::Int = 0
    match_ic::Float64 = 0.0
    normalized_ic::Float64 = 0.0
    core_ic::Float64 = 0.0
    score::Float64 = 0.0
    query_information::Float64 = 0.0
    search_information::Float64 = 0.0
end

@def_pprint mime_types="text/plain" base_show=true ComparisonResult

@enum _PositionKind::UInt8 begin
    _RESIDUE = 0
    _NTERMINUS = 1
    _CTERMINUS = 2
end

"""
    MatchFixMode

Fixed-position matching behavior used by CompariMotif:
- `MatchFixNone`: no fixed-position requirement.
- `MatchFixQueryFixed`: fixed query positions must have exact fixed matches.
- `MatchFixSearchFixed`: fixed search positions must have exact fixed matches.
- `MatchFixBothFixed`: enforce fixed-position matching on both motifs.
"""
@enum MatchFixMode::UInt8 begin
    MatchFixNone = 0
    MatchFixQueryFixed = 1
    MatchFixSearchFixed = 2
    MatchFixBothFixed = 3
end

struct _Position
    kind::_PositionKind
    mask::UInt32
end

struct _Token
    position::_Position
    min_repeat::Int
    max_repeat::Int
    canonical::String
end

struct _ParsedMotif
    original::String
    normalized::String
    tokens::Vector{_Token}
end

struct _MotifVariant
    positions::Vector{_Position}
    normalized::String
    information::Float64
end

"""
    ComparisonOptions

Reusable configuration object for CompariMotif comparisons.

Construct once with [`ComparisonOptions(; kwargs...)`](@ref) and reuse across
many `compare` calls.

# Keywords
- `alphabet::Symbol = :protein`: comparison alphabet (`:protein` or `:dna`).
- `min_shared_positions::Int = 2`: minimum number of matched, non-wildcard
  positions required for a hit.
- `normalized_ic_cutoff::Real = 0.5`: minimum normalized information content.
- `matchfix::Union{MatchFixMode, Symbol, AbstractString} = MatchFixNone`:
  fixed-position matching mode. Accepted symbol/string aliases are:
  `none`, `query_fixed` (`query`), `search_fixed` (`search`), `both_fixed` (`both`).
- `mismatches::Int = 0`: tolerated count of defined-position mismatches.
- `allow_ambiguous_overlap::Bool = true`: whether partial class overlaps are
  allowed as complex matches.
- `max_variants::Int = 10_000`: maximum expanded variants per motif.
"""
struct ComparisonOptions
    alphabet::Vector{Char}
    alphabet_index::Dict{Char, Int}
    alphabet_mask::UInt32
    log_base::Float64
    min_shared_positions::Int
    normalized_ic_cutoff::Float64
    matchfix::MatchFixMode
    mismatches::Int
    allow_ambiguous_overlap::Bool
    max_variants::Int
end

@def_pprint mime_types="text/plain" base_show=true ComparisonOptions

@enum _RelationshipType::UInt8 begin
    _REL_EXACT = 0
    _REL_VARIANT = 1
    _REL_DEGENERATE = 2
    _REL_COMPLEX = 3
end

@enum _RelationshipLength::UInt8 begin
    _LEN_MATCH = 0
    _LEN_PARENT = 1
    _LEN_SUBSEQUENCE = 2
    _LEN_OVERLAP = 3
end

struct _Candidate
    query_variant::_MotifVariant
    search_variant::_MotifVariant
    query_relationship_type::_RelationshipType
    query_relationship_length::_RelationshipLength
    search_relationship_type::_RelationshipType
    search_relationship_length::_RelationshipLength
    matched_pattern::String
    matched_positions::Int
    exact_fixed_matches::Int
    match_ic::Float64
    normalized_ic::Float64
    core_ic::Float64
    score::Float64
end

const _RELATIONSHIP_TYPE_WORDS = ("Exact", "Variant", "Degenerate", "Complex")
const _RELATIONSHIP_TYPE_CODES = ("e", "v", "d", "c")
const _RELATIONSHIP_LENGTH_WORDS = ("Match", "Parent", "Subsequence", "Overlap")
const _RELATIONSHIP_LENGTH_CODES = ("m", "p", "s", "o")

_coerce_matchfix(mode::MatchFixMode) = mode
function _coerce_matchfix(mode::Symbol)
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

    alphabet_chars = if alphabet === :protein
        _PROTEIN_ALPHABET
    elseif alphabet === :dna
        _DNA_ALPHABET
    else
        throw(ArgumentError("`alphabet` must be :protein or :dna."))
    end

    index = Dict{Char, Int}()
    for (i, aa) in enumerate(alphabet_chars)
        index[aa] = i
    end
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

_is_terminus(pos::_Position) = pos.kind !== _RESIDUE
function _is_wildcard(pos::_Position, options::ComparisonOptions)
    pos.kind === _RESIDUE && pos.mask == options.alphabet_mask
end
_is_fixed(pos::_Position) = pos.kind === _RESIDUE && count_ones(pos.mask) == 1

function _position_ic(pos::_Position, options::ComparisonOptions)
    if pos.kind !== _RESIDUE
        return 1.0
    end
    if _is_wildcard(pos, options)
        return 0.0
    end
    k = count_ones(pos.mask)
    return -log(k / length(options.alphabet)) / options.log_base
end

function _parse_repeat_quantifier(text::AbstractString, i::Int)
    if i > lastindex(text)
        return (1, 1, i)
    end
    opener = text[i]
    if opener != '{' && opener != '('
        return (1, 1, i)
    end
    closer = opener == '{' ? '}' : ')'
    close_idx = findnext(==(closer), text, nextind(text, i))
    close_idx === nothing &&
        throw(ArgumentError("Unclosed repeat quantifier in motif: $text"))
    raw = text[nextind(text, i):prevind(text, close_idx)]
    parts = split(raw, ',')
    repeat_min,
    repeat_max = if length(parts) == 1
        n = parse(Int, strip(parts[1]))
        (n, n)
    elseif length(parts) == 2
        n = parse(Int, strip(parts[1]))
        m = parse(Int, strip(parts[2]))
        (n, m)
    else
        throw(ArgumentError("Invalid repeat quantifier in motif: $text"))
    end
    if repeat_min < 0 || repeat_max < repeat_min
        throw(ArgumentError("Invalid repeat bounds in motif: $text"))
    end
    return (repeat_min, repeat_max, nextind(text, close_idx))
end

function _mask_from_char(char::Char, options::ComparisonOptions)
    aa = uppercase(char)
    idx = get(options.alphabet_index, aa, 0)
    idx == 0 && throw(ArgumentError("Unsupported residue '$char' for selected alphabet."))
    return UInt32(1) << (idx - 1)
end

function _class_mask(raw::AbstractString, options::ComparisonOptions)
    isempty(raw) && throw(ArgumentError("Empty character class is not allowed."))
    invert = startswith(raw, "^")
    body = invert ? raw[nextind(raw, firstindex(raw)):end] : raw
    isempty(body) && throw(ArgumentError("Empty negated character class is not allowed."))
    mask = UInt32(0)
    for char in body
        mask |= _mask_from_char(char, options)
    end
    if invert
        mask = options.alphabet_mask & ~mask
    end
    mask == 0 && throw(ArgumentError("Character class resolves to an empty set."))
    return mask
end

function _mask_to_chars(mask::UInt32, options::ComparisonOptions; as_lowercase::Bool = false)
    chars = Char[]
    for (i, aa) in enumerate(options.alphabet)
        if (mask & (UInt32(1) << (i - 1))) != 0
            push!(chars, as_lowercase ? Base.lowercase(aa) : aa)
        end
    end
    return chars
end

function _mask_to_symbol(mask::UInt32, options::ComparisonOptions;
        as_lowercase::Bool = false, wildcard_symbol::String = "x")
    if mask == options.alphabet_mask
        return as_lowercase ? Base.lowercase(wildcard_symbol) : wildcard_symbol
    end
    chars = _mask_to_chars(mask, options; as_lowercase = as_lowercase)
    length(chars) == 1 && return string(chars[1])
    return "[" * join(chars) * "]"
end

function _canonical_token(position::_Position, options::ComparisonOptions)
    if position.kind == _NTERMINUS
        return "^"
    elseif position.kind == _CTERMINUS
        return "\$"
    end
    return _mask_to_symbol(position.mask, options; as_lowercase = false, wildcard_symbol = "x")
end

function _parse_motif(motif::AbstractString, options::ComparisonOptions)
    stripped = strip(motif)
    isempty(stripped) && throw(ArgumentError("Motif cannot be empty."))

    tokens = _Token[]
    i = firstindex(stripped)
    while i <= lastindex(stripped)
        char = stripped[i]
        if isspace(char)
            i = nextind(stripped, i)
            continue
        end

        position = if char == '^'
            _Position(_NTERMINUS, 0)
        elseif char == '$'
            _Position(_CTERMINUS, 0)
        elseif char == 'x' || char == 'X' || char == '.'
            _Position(_RESIDUE, options.alphabet_mask)
        elseif char == '['
            close_idx = findnext(==(']'), stripped, nextind(stripped, i))
            close_idx === nothing &&
                throw(ArgumentError("Unclosed character class in motif: $motif"))
            class_raw = stripped[nextind(stripped, i):prevind(stripped, close_idx)]
            mask = _class_mask(class_raw, options)
            i = close_idx
            _Position(_RESIDUE, mask)
        else
            _Position(_RESIDUE, _mask_from_char(char, options))
        end

        repeat_min, repeat_max,
        next_i = _parse_repeat_quantifier(stripped, nextind(stripped, i))
        if position.kind !== _RESIDUE && (repeat_min != 1 || repeat_max != 1)
            throw(ArgumentError("Repeat quantifiers are not valid for termini in motif: $motif"))
        end
        if repeat_max == 0
            i = next_i
            continue
        end

        canonical = _canonical_token(position, options)
        if repeat_min == 1 && repeat_max == 1
            # keep canonical unchanged
        elseif repeat_min == repeat_max
            canonical *= "{" * string(repeat_min) * "}"
        else
            canonical *= "{" * string(repeat_min) * "," * string(repeat_max) * "}"
        end

        push!(tokens, _Token(position, repeat_min, repeat_max, canonical))
        i = next_i
    end

    isempty(tokens) &&
        throw(ArgumentError("Motif produced no positions after parsing: $motif"))
    normalized = join(getfield.(tokens, :canonical))
    return _ParsedMotif(String(motif), normalized, tokens)
end

function _variant_count(tokens::Vector{_Token})
    total = big(1)
    for token in tokens
        total *= (token.max_repeat - token.min_repeat + 1)
    end
    return total
end

function _expand_variants(parsed::_ParsedMotif, options::ComparisonOptions)
    nvariants = _variant_count(parsed.tokens)
    if nvariants > options.max_variants
        throw(ArgumentError("Motif $(parsed.original) expands to $nvariants variants, above max_variants=$(options.max_variants)."))
    end

    variants = _MotifVariant[]
    positions = _Position[]
    symbols = String[]

    function rec(ti::Int)
        if ti > length(parsed.tokens)
            isempty(positions) && return
            info = 0.0
            for pos in positions
                info += _position_ic(pos, options)
            end
            push!(variants, _MotifVariant(copy(positions), join(symbols), info))
            return
        end
        token = parsed.tokens[ti]
        token_symbol = _canonical_token(token.position, options)
        for repeat in token.min_repeat:token.max_repeat
            append_count = 0
            for _ in 1:repeat
                push!(positions, token.position)
                push!(symbols, token_symbol)
                append_count += 1
            end
            rec(ti + 1)
            for _ in 1:append_count
                pop!(positions)
                pop!(symbols)
            end
        end
    end
    rec(1)
    return variants
end

function _match_symbol(
        qpos::_Position,
        spos::_Position,
        intersection::UInt32,
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
        return "x"
    end
    if mismatch
        union_mask = qpos.mask | spos.mask
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
    union_mask = qpos.mask | spos.mask
    return _mask_to_symbol(union_mask, options; as_lowercase = true, wildcard_symbol = "x")
end

function _relationship_type_from_flags(has_variant::Bool, has_degenerate::Bool, has_complex::Bool)
    if has_complex || (has_variant && has_degenerate)
        return _REL_COMPLEX
    elseif has_variant
        return _REL_VARIANT
    elseif has_degenerate
        return _REL_DEGENERATE
    else
        return _REL_EXACT
    end
end

function _reverse_type(relationship_type::_RelationshipType)
    if relationship_type == _REL_VARIANT
        return _REL_DEGENERATE
    elseif relationship_type == _REL_DEGENERATE
        return _REL_VARIANT
    end
    return relationship_type
end

function _relationship_length(qlen::Int, slen::Int, overlap::Int)
    if overlap == qlen && qlen == slen
        return _LEN_MATCH
    elseif overlap == slen && qlen > slen
        return _LEN_PARENT
    elseif overlap == qlen && slen > qlen
        return _LEN_SUBSEQUENCE
    else
        return _LEN_OVERLAP
    end
end

function _reverse_length(length_type::_RelationshipLength)
    if length_type == _LEN_PARENT
        return _LEN_SUBSEQUENCE
    elseif length_type == _LEN_SUBSEQUENCE
        return _LEN_PARENT
    end
    return length_type
end

function _relationship_word(relationship_type::_RelationshipType, length_type::_RelationshipLength)
    return _RELATIONSHIP_TYPE_WORDS[Int(relationship_type) + 1] * " " *
           _RELATIONSHIP_LENGTH_WORDS[Int(length_type) + 1]
end

function _relationship_code(relationship_type::_RelationshipType, length_type::_RelationshipLength)
    return _RELATIONSHIP_TYPE_CODES[Int(relationship_type) + 1] * "-" *
           _RELATIONSHIP_LENGTH_CODES[Int(length_type) + 1]
end

function _compare_positions(qpos::_Position, spos::_Position, options::ComparisonOptions)
    if qpos.kind !== _RESIDUE || spos.kind !== _RESIDUE
        if qpos.kind == spos.kind
            return (hard_mismatch = false, mismatch = false,
                relation = _REL_EXACT, intersection = UInt32(0), ic = 1.0,
                contributes_position = true, exact_fixed = false)
        end
        return (hard_mismatch = true, mismatch = true,
            relation = _REL_COMPLEX, intersection = UInt32(0), ic = 0.0,
            contributes_position = false, exact_fixed = false)
    end

    intersection = qpos.mask & spos.mask
    if intersection == 0
        return (hard_mismatch = false, mismatch = true,
            relation = _REL_COMPLEX, intersection = intersection, ic = 0.0,
            contributes_position = false, exact_fixed = false)
    end

    relation = if qpos.mask == spos.mask
        _REL_EXACT
    elseif (qpos.mask & ~spos.mask) == 0
        _REL_VARIANT
    elseif (spos.mask & ~qpos.mask) == 0
        _REL_DEGENERATE
    elseif options.allow_ambiguous_overlap
        _REL_COMPLEX
    else
        return (hard_mismatch = true, mismatch = true, relation = _REL_COMPLEX,
            intersection = intersection, ic = 0.0,
            contributes_position = false, exact_fixed = false)
    end

    ic = min(_position_ic(qpos, options), _position_ic(spos, options))
    contributes = !_is_wildcard(qpos, options) && !_is_wildcard(spos, options)
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

function _query_fixed_required(mode::MatchFixMode)
    mode == MatchFixQueryFixed || mode == MatchFixBothFixed
end
function _search_fixed_required(mode::MatchFixMode)
    mode == MatchFixSearchFixed || mode == MatchFixBothFixed
end

function _evaluate_alignment(
        query_variant::_MotifVariant,
        search_variant::_MotifVariant,
        shift::Int,
        options::ComparisonOptions
)
    qlen = length(query_variant.positions)
    slen = length(search_variant.positions)
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
                return nothing
            end
        end
        if _search_fixed_required(options.matchfix)
            if _is_fixed(spos) && !cmp.exact_fixed
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

function _compare_parsed(parsed_query::_ParsedMotif, parsed_search::_ParsedMotif, options::ComparisonOptions)
    query_variants = _expand_variants(parsed_query, options)
    search_variants = _expand_variants(parsed_search, options)
    best = nothing

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
        query_relationship_code = _relationship_code(
            best.query_relationship_type, best.query_relationship_length),
        search_relationship_code = _relationship_code(
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

"""
    normalize_motif(motif::AbstractString; alphabet::Symbol = :protein) -> String

Parse and canonicalize a motif expression into a deterministic representation.
Supported syntax includes fixed residues, bracket classes (including negation),
`x`/`.` wildcards, `^`/`\$` termini, and `{m,n}` (or `(m,n)`) repeat quantifiers.
"""
function normalize_motif(motif::AbstractString; alphabet::Symbol = :protein)
    options = ComparisonOptions(; alphabet, min_shared_positions = 1, normalized_ic_cutoff = 0.0)
    return _parse_motif(motif, options).normalized
end

"""
    compare(a::AbstractString, b::AbstractString, options::ComparisonOptions) -> ComparisonResult

Compare two motifs and return the best relationship according to the
CompariMotif scoring scheme described in Edwards et al. (2008).
"""
function compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)
    parsed_a = _parse_motif(a, options)
    parsed_b = _parse_motif(b, options)
    return _compare_parsed(parsed_a, parsed_b, options)
end

"""
    compare(motifs::AbstractVector{<:AbstractString},
            db::AbstractVector{<:AbstractString},
            options::ComparisonOptions) -> Matrix{ComparisonResult}

Compute all pairwise comparisons between query motifs and database motifs.
The result matrix has size `(length(motifs), length(db))`.
"""
function compare(
        motifs::AbstractVector{<:AbstractString},
        db::AbstractVector{<:AbstractString},
        options::ComparisonOptions
)
    parsed_queries = [_parse_motif(motif, options) for motif in motifs]
    parsed_db = [_parse_motif(motif, options) for motif in db]

    results = Matrix{ComparisonResult}(undef, length(parsed_queries), length(parsed_db))
    for i in eachindex(parsed_queries)
        for j in eachindex(parsed_db)
            results[i, j] = _compare_parsed(parsed_queries[i], parsed_db[j], options)
        end
    end
    return results
end

"""
    compare(motifs::AbstractVector{<:AbstractString},
            options::ComparisonOptions) -> Matrix{ComparisonResult}

Convenience method for all-vs-all motif comparison.
"""
function compare(motifs::AbstractVector{<:AbstractString}, options::ComparisonOptions)
    compare(motifs, motifs, options)
end

"""
    write_results_tsv(path, motifs, db, results)

Write pairwise comparison results to a deterministic TSV file.
"""
function write_results_tsv(path::AbstractString,
        motifs::AbstractVector{<:AbstractString},
        db::AbstractVector{<:AbstractString},
        results::Matrix{ComparisonResult})
    size(results) == (length(motifs), length(db)) ||
        throw(ArgumentError("`results` dimensions must match `(length(motifs), length(db))`."))

    open(path, "w") do io
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

end
