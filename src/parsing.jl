"""
    _is_terminus(pos::_Position)::Bool

Return `true` when `pos` is a terminus anchor (`^` or `\$`).
"""
_is_terminus(pos::_Position) = pos.kind !== _RESIDUE

"""
    _is_wildcard(pos::_Position, wildcard_mask::ResidueMask)::Bool

Return `true` when `pos` matches all residues in the selected alphabet.
"""
function _is_wildcard(pos::_Position, wildcard_mask::ResidueMask)
    pos.kind === _RESIDUE && pos.mask == wildcard_mask
end

"""
    _is_fixed(pos::_Position)::Bool

Return `true` when `pos` encodes exactly one residue.
"""
_is_fixed(pos::_Position) = pos.kind === _RESIDUE && is_fixed(ResidueClass(pos.mask))

"""
    _position_ic(pos::_Position, spec::_AlphabetSpec)::Float64

Compute information content for one parsed position.
"""
function _position_ic(pos::_Position, spec::_AlphabetSpec)
    return _position_ic(pos, spec, _uniform_residue_frequency_vector(spec))
end

function _position_ic(
        pos::_Position,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    # Termini contribute as defined anchors.
    if pos.kind !== _RESIDUE
        return 1.0
    end
    # Wildcards carry no residue specificity and therefore no IC contribution.
    if _is_wildcard(pos, spec.mask)
        return 0.0
    end
    # Information depends on the total background mass of the residue set.
    mass = 0.0
    for i in eachindex(residue_frequencies)
        if _mask_has_index(pos.mask, i)
            mass += residue_frequencies[i]
        end
    end
    return -log(mass) / spec.log_base
end

"""
    _parse_repeat_quantifier(text::AbstractString, i::Int)::Tuple{Int, Int, Int}

Parse optional repeat quantifier at index `i`, returning `(min, max, next_index)`.
"""
function _parse_repeat_quantifier(text::AbstractString, i::Int)
    # No suffix quantifier -> implicit {1}.
    if i > lastindex(text)
        return (1, 1, i)
    end
    opener = text[i]
    if opener != '{'
        return (1, 1, i)
    end
    close_idx = findnext(==('}'), text, nextind(text, i))
    close_idx === nothing &&
        throw(ArgumentError("Unclosed repeat quantifier in motif: $text"))
    raw = text[nextind(text, i):prevind(text, close_idx)]
    parts = split(raw, ',')
    # Supported forms:
    # - "{n}" -> exact repeat
    # - "{m,n}" -> range repeat
    repeat_min,
    repeat_max = if length(parts) == 1
        n = parse(Int, strip(parts[1]))
        (max(n, 1), max(n, 1))
    elseif length(parts) == 2
        n = parse(Int, strip(parts[1]))
        m = parse(Int, strip(parts[2]))
        (n, m)
    else
        throw(ArgumentError("Invalid repeat quantifier in motif: $text"))
    end
    if repeat_min < 0 || repeat_max < 0
        throw(ArgumentError("Invalid repeat bounds in motif: $text"))
    end
    return (repeat_min, repeat_max, nextind(text, close_idx))
end

# -----------------------------------------------------------------------------
# Residue masks
# -----------------------------------------------------------------------------
# Internal residue sets are represented as bit masks.
# - A single residue sets one bit.
# - Character classes set multiple bits.
# - Wildcards use the full alphabet mask.
# - Negated classes are computed as complement within the selected alphabet.
#
# This gives stable set operations for exact/subset/superset/overlap checks.

"""
    _mask_from_char(char::Char, spec::_AlphabetSpec)::ResidueMask

Return the residue mask for one alphabet character.
"""
function _mask_from_char(char::Char, spec::_AlphabetSpec)
    # Matching is case-insensitive at parse time.
    aa = uppercase(char)
    idx = get(spec.index, aa, 0)
    idx == 0 && throw(ArgumentError("Unsupported residue '$char' for selected alphabet."))
    # The first residue in alphabet uses bit 0, second uses bit 1, etc.
    return _residue_mask_bit(idx)
end

"""
    _class_mask(raw::AbstractString, spec::_AlphabetSpec)::ResidueMask

Parse a bracket class body into a residue mask.
"""
function _class_mask(raw::AbstractString, spec::_AlphabetSpec)
    isempty(raw) && throw(ArgumentError("Empty character class is not allowed."))
    # `[^...]` syntax means "all residues except listed residues".
    invert = startswith(raw, "^")
    body = invert ? raw[nextind(raw, firstindex(raw)):end] : raw
    isempty(body) && throw(ArgumentError("Empty negated character class is not allowed."))
    mask = ResidueMask(0)
    for char in body
        # Class set is the union of member residue bits.
        mask |= _mask_from_char(char, spec)
    end
    if invert
        # Complement inside alphabet domain only.
        mask = _mask_complement(mask, spec)
    end
    mask == 0 && throw(ArgumentError("Character class resolves to an empty set."))
    return mask
end

"""
    _mask_to_chars(mask::ResidueMask, spec::_AlphabetSpec; as_lowercase = false)::Vector{Char}

Materialize residues represented by a mask in canonical alphabet order.
"""
function _mask_to_chars(mask::ResidueMask, spec::_AlphabetSpec; as_lowercase::Bool = false)
    chars = Char[]
    for (i, aa) in enumerate(spec.chars)
        # Emit residues in canonical alphabet order for deterministic normalization.
        if _mask_has_index(mask, i)
            push!(chars, as_lowercase ? Base.lowercase(aa) : aa)
        end
    end
    return chars
end

"""
    _mask_to_symbol(mask::ResidueMask, spec::_AlphabetSpec; as_lowercase = false, wildcard_symbol = ".")::String

Render one residue mask as canonical motif syntax.
"""
function _mask_to_symbol(mask::ResidueMask, spec::_AlphabetSpec;
        as_lowercase::Bool = false, wildcard_symbol::String = ".")
    # Full mask is represented as wildcard, not as long explicit class.
    if mask == spec.mask
        return as_lowercase ? Base.lowercase(wildcard_symbol) : wildcard_symbol
    end
    chars = _mask_to_chars(mask, spec; as_lowercase = as_lowercase)
    length(chars) == 1 && return string(chars[1])
    return "[" * join(chars) * "]"
end

"""
    _canonical_token(position::_Position, spec::_AlphabetSpec)::String

Render one parsed position into deterministic canonical motif syntax.
"""
function _canonical_token(position::_Position, spec::_AlphabetSpec)
    # Keep termini canonicalized exactly as anchors.
    if position.kind == _NTERMINUS
        return "^"
    elseif position.kind == _CTERMINUS
        return "\$"
    end
    # Residues/classes canonicalize through mask representation.
    return _mask_to_symbol(position.mask, spec)
end

# -----------------------------------------------------------------------------
# Motif parser
# -----------------------------------------------------------------------------
# The parser keeps syntax handling deterministic:
# - fixed residues, wildcard tokens (`x`, `X`, `.`), classes, negated classes,
#   and termini are parsed into `_Position` tokens;
# - optional quantifiers `{n}` and `{m,n}` are attached to the preceding token;
# - grouping and alternation (`(...)`, `|`) are expanded into independent branch
#   token lists;
# - leading/trailing whitespace is trimmed, while the first internal whitespace
#   ends parsing to mirror the upstream oracle;
# - canonical text is rebuilt from masks and quantifiers to normalize input.

"""
    _normalized_from_tokens(tokens::Vector{_Token})::String

Join canonical token renderings into one normalized motif string.

# Examples
```jldoctest
julia> using CompariMotif

julia> options = CompariMotif.ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> parsed = CompariMotif._parse_motif("r[kR].{0,1}l", options);

julia> CompariMotif._normalized_from_tokens(parsed.tokens)
"R[RK].{0,1}L"
```
"""
function _normalized_from_tokens(tokens::Vector{_Token})
    return join(getfield.(tokens, :canonical))
end

"""
    _variant_limit_error(motif::AbstractString, nvariants::BigInt, max_variants::Int)::ArgumentError

Construct the deterministic error used when grouped alternation or repeat
expansion would exceed `max_variants`.
"""
function _variant_limit_error(motif::AbstractString, nvariants::BigInt, max_variants::Int)
    return ArgumentError(
        "Motif $motif expands to $nvariants variants, above max_variants=$max_variants."
    )
end

"""
    _oracle_parse_window(motif::AbstractString)::String

Trim leading/trailing whitespace, then stop at the first internal whitespace to
mirror the upstream parser's effective input window.

# Examples
```jldoctest
julia> using CompariMotif

julia> CompariMotif._oracle_parse_window("  A(K|Q) LI  ")
"A(K|Q)"
```
"""
function _oracle_parse_window(motif::AbstractString)
    stripped = strip(motif)
    isempty(stripped) && return stripped

    i = firstindex(stripped)
    while i <= lastindex(stripped)
        if isspace(stripped[i])
            return stripped[firstindex(stripped):prevind(stripped, i)]
        end
        i = nextind(stripped, i)
    end
    return stripped
end

"""
    _combine_concat(lhs, rhs, motif, max_variants)::Vector{String}

Form the Cartesian concatenation of two alternative sets while enforcing the
variant expansion limit for the original `motif`.

# Examples
```jldoctest
julia> using CompariMotif

julia> CompariMotif._combine_concat(["A", "B"], ["C", "D"], "demo", 10)
4-element Vector{String}:
 "AC"
 "AD"
 "BC"
 "BD"
```
"""
function _combine_concat(
        lhs::Vector{String},
        rhs::Vector{String},
        motif::AbstractString,
        max_variants::Int
)
    total = big(length(lhs)) * big(length(rhs))
    total > max_variants && throw(_variant_limit_error(motif, total, max_variants))

    out = String[]
    sizehint!(out, Int(total))
    for left in lhs, right in rhs

        push!(out, left * right)
    end
    return out
end

"""
    _parse_expr_alternatives(text::AbstractString, i::Int, max_variants::Int)

Parse the expression starting at `i` into all top-level alternation branches,
stopping at `')'` or the end of `text`.
"""
function _parse_expr_alternatives(text::AbstractString, i::Int, max_variants::Int)
    terms = String[]
    seq = [""]

    while true
        if i > lastindex(text) || text[i] == ')'
            append!(terms, seq)
            return terms, i
        end
        if text[i] == '|'
            isempty(seq) && throw(ArgumentError("Empty alternation branch in motif: $text"))
            append!(terms, seq)
            seq = [""]
            i = nextind(text, i)
            continue
        end

        atom_alts, i = _parse_atom_alternatives(text, i, max_variants)
        isempty(atom_alts) && return (String[], i)
        seq = _combine_concat(seq, atom_alts, text, max_variants)
    end
end

"""
    _extract_quantifier(text::AbstractString, i::Int)::Tuple{String, Int}

Return the raw `{...}` quantifier starting at `i`, or `("", i)` when no repeat
suffix is present.
"""
function _extract_quantifier(text::AbstractString, i::Int)
    if i > lastindex(text) || text[i] != '{'
        return "", i
    end
    _, _, next_i = _parse_repeat_quantifier(text, i)
    return String(text[i:prevind(text, next_i)]), next_i
end

"""
    _parse_atom_alternatives(text::AbstractString, i::Int, max_variants::Int)

Parse one atom at `i`, including bracket classes, grouped alternation, and any
trailing repeat quantifier, then return its concrete textual alternatives.
"""
function _parse_atom_alternatives(text::AbstractString, i::Int, max_variants::Int)
    char = text[i]

    if char == '['
        close_idx = findnext(==(']'), text, nextind(text, i))
        close_idx === nothing &&
            throw(ArgumentError("Unclosed character class in motif: $text"))
        atom = String(text[i:close_idx])
        quant, next_i = _extract_quantifier(text, nextind(text, close_idx))
        return [atom * quant], next_i
    end
    if char == '('
        inner_alts, next_i = _parse_expr_alternatives(text, nextind(text, i), max_variants)
        isempty(inner_alts) && return (String[], next_i)
        next_i > lastindex(text) &&
            throw(ArgumentError("Unclosed grouping parenthesis in motif: $text"))
        text[next_i] == ')' || throw(ArgumentError("Malformed grouping in motif: $text"))
        quant, tail_i = _extract_quantifier(text, nextind(text, next_i))
        if isempty(quant)
            return inner_alts, tail_i
        end
        return [alt * quant for alt in inner_alts], tail_i
    end
    quant, next_i = _extract_quantifier(text, nextind(text, i))
    return [string(char) * quant], next_i
end

"""
    _expand_grouped_motif(text::AbstractString, max_variants::Int)::Vector{String}

Expand grouped alternation syntax into canonical branch strings before token
parsing.

# Examples
```jldoctest
julia> using CompariMotif

julia> CompariMotif._expand_grouped_motif("A(K|Q)L", 10)
2-element Vector{String}:
 "AKL"
 "AQL"
```
"""
function _expand_grouped_motif(text::AbstractString, max_variants::Int)
    alts, next_i = _parse_expr_alternatives(text, firstindex(text), max_variants)
    isempty(alts) && return alts
    next_i <= lastindex(text) && text[next_i] != ')' &&
        throw(ArgumentError("Unexpected trailing content in motif: $text"))

    cleaned = String[]
    for alt in alts
        stripped = strip(alt)
        isempty(stripped) &&
            throw(ArgumentError("Empty alternation branch in motif: $text"))
        push!(cleaned, stripped)
    end
    return cleaned
end

function _parse_linear_tokens(motif::AbstractString, spec::_AlphabetSpec)
    tokens = _Token[]
    i = firstindex(motif)
    while i <= lastindex(motif)
        char = motif[i]

        # Parse a single symbolic unit (position/anchor/class/wildcard).
        position = if char == '^'
            _Position(_NTERMINUS, 0)
        elseif char == '$'
            _Position(_CTERMINUS, 0)
        elseif char == 'x' || char == 'X' || char == '.'
            # Wildcards match any residue in the selected alphabet.
            _Position(_RESIDUE, spec.mask)
        elseif char == '['
            # Class token includes everything until the first closing bracket.
            close_idx = findnext(==(']'), motif, nextind(motif, i))
            close_idx === nothing &&
                throw(ArgumentError("Unclosed character class in motif: $motif"))
            class_raw = motif[nextind(motif, i):prevind(motif, close_idx)]
            mask = _class_mask(class_raw, spec)
            i = close_idx
            _Position(_RESIDUE, mask)
        else
            _Position(_RESIDUE, _mask_from_char(char, spec))
        end

        repeat_min, repeat_max,
        next_i = _parse_repeat_quantifier(motif, nextind(motif, i))

        canonical = _canonical_token(position, spec)
        if repeat_min == 1 && repeat_max == 1
            # no quantifier suffix
        elseif repeat_min == repeat_max
            canonical *= "{" * string(repeat_min) * "}"
        else
            canonical *= "{" * string(repeat_min) * "," * string(repeat_max) * "}"
        end

        push!(tokens, _Token(position, repeat_min, repeat_max, canonical))
        # Advance to the next symbol after an optional quantifier.
        i = next_i
    end

    isempty(tokens) &&
        throw(ArgumentError("Motif produced no positions after parsing: $motif"))
    return tokens
end

"""
    _parse_motif(motif::AbstractString, options::ComparisonOptions)::_ParsedMotif

Parse one motif string into canonical internal representation.
"""
function _parse_motif(motif::AbstractString, options::ComparisonOptions)
    parse_window = _oracle_parse_window(motif)
    isempty(parse_window) && throw(ArgumentError("Motif cannot be empty."))

    spec = _alphabet_spec(options.alphabet)
    branches = _expand_grouped_motif(parse_window, options.max_variants)
    isempty(branches) &&
        return _ParsedMotif(String(motif), parse_window, _Token[], Vector{_Token}[])
    alternatives = Vector{_Token}[]
    for branch in branches
        branch_tokens = _parse_linear_tokens(branch, spec)
        push!(alternatives, branch_tokens)
    end
    if length(alternatives) > 1
        _canonicalize_alternatives!(alternatives, spec)
    end
    normalized_branches = [_normalized_from_tokens(branch_tokens)
                           for branch_tokens in alternatives]
    normalized = if length(normalized_branches) == 1
        first(normalized_branches)
    else
        join(["(" * branch * ")" for branch in normalized_branches], "|")
    end
    return _ParsedMotif(String(motif), normalized, first(alternatives), alternatives)
end

"""
    _variant_count(tokens::Vector{_Token})::BigInt

Return the number of expanded variants implied by repeat ranges.
"""
function _variant_count(tokens::Vector{_Token})
    # Product of each token's repeat choice count, computed in BigInt to avoid overflow.
    total = big(1)
    for token in tokens
        total *= (token.max_repeat - token.min_repeat + 1)
    end
    return total
end

"""
    _expand_variants(parsed::_ParsedMotif, options::ComparisonOptions, spec::_AlphabetSpec)::Vector{_MotifVariant}

Expand ranged-repeat motifs into concrete variant sequences.
"""
function _expand_variants(parsed::_ParsedMotif, options::ComparisonOptions, spec::_AlphabetSpec)
    return _expand_variants(parsed, options, spec, _residue_frequency_vector(options, spec))
end

function _expand_variants(
        parsed::_ParsedMotif,
        options::ComparisonOptions,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    nvariants = big(0)
    for tokens in parsed.alternatives
        nvariants += _variant_count(tokens)
    end
    if nvariants > options.max_variants
        throw(_variant_limit_error(parsed.original, nvariants, options.max_variants))
    end

    variants = _MotifVariant[]

    for tokens in parsed.alternatives
        append!(variants, _expand_branch_variants(tokens, spec, residue_frequencies))
    end
    return variants
end

function _emit_branch_variant!(
        variants::Vector{_MotifVariant},
        tokens::Vector{_Token},
        repeat_counts::Vector{Int},
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    # Build one concrete variant from the currently chosen repeat count per token.
    positions = _Position[]
    symbols = String[]
    for (token, repeat) in zip(tokens, repeat_counts)
        # Reuse the token's canonical symbol so variant text stays deterministic.
        token_symbol = _canonical_token(token.position, spec)
        for _ in 1:repeat
            push!(positions, token.position)
            push!(symbols, token_symbol)
        end
    end
    # If all repeats are zero, this branch contributes no concrete variant.
    isempty(positions) && return nothing

    # Variant information is the sum of per-position information content.
    info = 0.0
    for pos in positions
        info += _position_ic(pos, spec, residue_frequencies)
    end
    push!(variants, _MotifVariant(positions, join(symbols), info))
    return nothing
end

function _enumerate_branch_repeat_counts!(
        variants::Vector{_MotifVariant},
        tokens::Vector{_Token},
        repeat_counts::Vector{Int},
        ti::Int,
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    # Base case: every token has a chosen repeat count, so emit one variant.
    ti > length(tokens) &&
        return _emit_branch_variant!(
            variants, tokens, repeat_counts, spec, residue_frequencies)
    token = tokens[ti]
    # Enumerate this token's allowed repeat counts, then recurse to later tokens.
    # This preserves oracle encounter order: earlier tokens vary more slowly,
    # later tokens vary faster.
    for repeat in token.min_repeat:token.max_repeat
        repeat_counts[ti] = repeat
        _enumerate_branch_repeat_counts!(
            variants,
            tokens,
            repeat_counts,
            ti + 1,
            spec,
            residue_frequencies
        )
    end
    return nothing
end

function _expand_branch_variants(
        tokens::Vector{_Token},
        spec::_AlphabetSpec,
        residue_frequencies::AbstractVector{<:Real}
)
    # Collect concrete variants for one alternation branch.
    variants = _MotifVariant[]
    # Mutable traversal state: selected repeat count for each token index.
    repeat_counts = zeros(Int, length(tokens))
    _enumerate_branch_repeat_counts!(
        variants,
        tokens,
        repeat_counts,
        1,
        spec,
        residue_frequencies
    )
    return variants
end

"""
    _canonicalize_alternatives!(alternatives::Vector{Vector{_Token}}, spec::_AlphabetSpec)

Sort grouped branch alternatives in place using the deterministic
oracle-backed branch-order rules.
"""
function _canonicalize_alternatives!(
        alternatives::Vector{Vector{_Token}},
        spec::_AlphabetSpec
)
    length(alternatives) <= 1 && return alternatives

    # Decorate each branch with temporary metadata
    # `(branch_index, tokens, order_tokens)`, sort using `order_tokens`,
    # then undecorate by writing only `tokens` back into `alternatives`.
    # This avoids recomputing order keys on every pairwise comparison.
    branch_entries = NamedTuple[]
    for (branch_index, tokens) in enumerate(alternatives)
        push!(branch_entries,
            (;
                branch_index,
                tokens,
                order_tokens = _branch_order_tokens(tokens)
            ))
    end

    sort!(branch_entries; lt = (left, right) -> _branch_tokens_isless(left, right, spec))
    empty!(alternatives)
    append!(alternatives, getfield.(branch_entries, :tokens))
    return alternatives
end

"""
    _branch_order_tokens(tokens::Vector{_Token})

Expand branch tokens into a comparison stream that separates mandatory copies
from optional repeat tails.
"""
function _branch_order_tokens(tokens::Vector{_Token})
    order_tokens = NamedTuple{(:position, :min_repeat, :max_repeat)}[]
    for token in tokens
        # `(min_repeat, max_repeat)` encode repeat bounds for this emitted unit.
        # Emit one exact `(1,1)` entry per mandatory copy.
        # Example: `R{2}` contributes two entries for `R`, both `(min_repeat=1, max_repeat=1)`.
        for _ in 1:token.min_repeat
            push!(order_tokens, (position = token.position, min_repeat = 1, max_repeat = 1))
        end
        # Any extra repeat capacity is stored as one optional tail `(0, extra)`.
        # Here `min_repeat=0` means the optional tail may contribute zero extra copies.
        # Example: `R{2,4}` contributes two exact entries plus one optional `(0,2)` tail.
        # Example: `R{0,3}` contributes only one optional `(0,3)` tail.
        extra_repeat = token.max_repeat - token.min_repeat
        if extra_repeat > 0
            push!(order_tokens,
                (position = token.position, min_repeat = 0, max_repeat = extra_repeat))
        end
    end
    return order_tokens
end

"""
    _position_order_category(pos::_Position, spec::_AlphabetSpec)

Return the coarse sort category used for branch ordering:
fixed residue < wildcard < residue class < N-terminus < C-terminus.
Returned rank values are:
- 0 fixed residue
- 1 wildcard
- 2 residue class
- 3 N-terminus
- 4 C-terminus
"""
@inline function _position_order_category(pos::_Position, spec::_AlphabetSpec)
    pos.kind === _NTERMINUS && return 3
    pos.kind === _CTERMINUS && return 4
    _is_fixed(pos) && return 0
    _is_wildcard(pos, spec.mask) && return 1
    return 2
end

"""
    _position_order_text(pos::_Position, spec::_AlphabetSpec)

Return the textual tie-break key for residue positions; anchors use `""`.
"""
@inline function _position_order_text(pos::_Position, spec::_AlphabetSpec)
    pos.kind === _RESIDUE || return ""
    return String(_mask_to_chars(pos.mask, spec))
end

"""
    _order_token_repeat_key(token)::Tuple{Int, Int}

Return the sort key `(repeat_kind_rank, max_repeat)` for one order token.

`repeat_kind_rank` values:
- `0`: exact repeat (`min_repeat == max_repeat`), e.g. `{3}`.
- `1`: ranged repeat (`min_repeat != max_repeat`), e.g. `{2,4}` or optional tail `(0,2)`.
"""
@inline function _order_token_repeat_key(token)::Tuple{Int, Int}
    is_exact = token.min_repeat == token.max_repeat
    return (is_exact ? 0 : 1, token.max_repeat)
end

"""
    _same_position(left::_Position, right::_Position)

Return `true` when two parsed positions encode the same anchor/residue mask.
"""
@inline function _same_position(left::_Position, right::_Position)
    return left.kind === right.kind && left.mask == right.mask
end

"""
    _is_exact_order_token(token)::Bool

Return `true` for the normalized exact-repeat unit `(1, 1)`.
"""
@inline function _is_exact_order_token(token)::Bool
    return token.min_repeat == 1 && token.max_repeat == 1
end

"""
    _is_optional_repeat_order_token(token)::Bool

Return `true` for normalized optional repeat tails `(0, n)` where `n > 0`.
"""
@inline function _is_optional_repeat_order_token(token)::Bool
    return token.min_repeat == 0 && token.max_repeat > 0
end

"""
    _repeat_overlap_ordering(left_order_tokens, right_order_tokens, idx::Int)

Apply the special tie-break for adjacent repeated positions that differ only by
`exact` versus `optional tail` representation at `idx`.
"""
function _repeat_overlap_ordering(
        left_order_tokens,
        right_order_tokens,
        idx::Int
)
    idx == 1 && return nothing

    left_token = left_order_tokens[idx]
    right_token = right_order_tokens[idx]
    _same_position(left_token.position, right_token.position) || return nothing

    left_exact = _is_exact_order_token(left_token)
    right_exact = _is_exact_order_token(right_token)
    left_optional = _is_optional_repeat_order_token(left_token)
    right_optional = _is_optional_repeat_order_token(right_token)
    if !(left_optional && right_exact) && !(left_exact && right_optional)
        return nothing
    end

    prev_left = left_order_tokens[idx - 1]
    prev_right = right_order_tokens[idx - 1]
    if !_same_position(prev_left.position, left_token.position) ||
       !_same_position(prev_right.position, right_token.position)
        return nothing
    end

    # Oracle-compatible preference: when one side has an exact extra copy and the
    # other side encodes that same copy as optional tail, keep the shorter side first.
    return left_optional && right_exact
end

"""
    _order_token_isless(left, right, spec::_AlphabetSpec)

Return `true` when `left` should sort before `right` by position category,
position text, then repeat key.
"""
function _order_token_isless(left, right, spec::_AlphabetSpec)
    # Primary key: broad position type.
    left_category = _position_order_category(left.position, spec)
    right_category = _position_order_category(right.position, spec)
    if left_category != right_category
        return left_category < right_category
    end

    # Secondary key: canonical residue text for deterministic class ordering.
    left_text = _position_order_text(left.position, spec)
    right_text = _position_order_text(right.position, spec)
    if left_text != right_text
        return left_text < right_text
    end

    # Tertiary key: exact repeats before ranged repeats.
    left_repeat = _order_token_repeat_key(left)
    right_repeat = _order_token_repeat_key(right)
    if left_repeat != right_repeat
        return left_repeat < right_repeat
    end
    return false
end

"""
    _branch_token_sequence_isless(left_order_tokens, right_order_tokens, spec::_AlphabetSpec)

Lexicographically compare two expanded branch token streams.
"""
function _branch_token_sequence_isless(left_order_tokens, right_order_tokens, spec::_AlphabetSpec)
    nshared = min(length(left_order_tokens), length(right_order_tokens))
    for idx in 1:nshared
        left_token = left_order_tokens[idx]
        right_token = right_order_tokens[idx]
        # First apply the repeat-overlap exception before normal token ordering.
        repeat_overlap_order = _repeat_overlap_ordering(left_order_tokens, right_order_tokens, idx)
        repeat_overlap_order === nothing || return repeat_overlap_order
        if _order_token_isless(left_token, right_token, spec)
            return true
        end
        if _order_token_isless(right_token, left_token, spec)
            return false
        end
    end
    return length(left_order_tokens) < length(right_order_tokens)
end

"""
    _branch_tokens_isless(left, right, spec::_AlphabetSpec)

Compare two decorated branch entries, falling back to original branch index to
keep sorting stable when order keys tie.
"""
function _branch_tokens_isless(left, right, spec::_AlphabetSpec)
    if _branch_token_sequence_isless(left.order_tokens, right.order_tokens, spec)
        return true
    end
    if _branch_token_sequence_isless(right.order_tokens, left.order_tokens, spec)
        return false
    end
    return left.branch_index < right.branch_index
end
