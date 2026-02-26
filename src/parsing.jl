# Utility predicates reused by parser and alignment.
_is_terminus(pos::_Position) = pos.kind !== _RESIDUE
function _is_wildcard(pos::_Position, options::ComparisonOptions)
    pos.kind === _RESIDUE && pos.mask == options.alphabet_mask
end
_is_fixed(pos::_Position) = pos.kind === _RESIDUE && count_ones(pos.mask) == 1

function _position_ic(pos::_Position, options::ComparisonOptions)
    # Termini contribute as defined anchors.
    if pos.kind !== _RESIDUE
        return 1.0
    end
    # Wildcards carry no residue specificity and therefore no IC contribution.
    if _is_wildcard(pos, options)
        return 0.0
    end
    # For a residue set of size k in alphabet size N, use log_N(1/(k/N)).
    k = count_ones(pos.mask)
    return -log(k / length(options.alphabet)) / options.log_base
end

function _parse_repeat_quantifier(text::AbstractString, i::Int)
    # No suffix quantifier -> implicit {1}.
    if i > lastindex(text)
        return (1, 1, i)
    end
    opener = text[i]
    # Accept both brace and parenthesis forms for compatibility.
    if opener != '{' && opener != '('
        return (1, 1, i)
    end
    closer = opener == '{' ? '}' : ')'
    close_idx = findnext(==(closer), text, nextind(text, i))
    close_idx === nothing &&
        throw(ArgumentError("Unclosed repeat quantifier in motif: $text"))
    raw = text[nextind(text, i):prevind(text, close_idx)]
    parts = split(raw, ',')
    # Supported forms:
    # - "{n}" or "(n)" -> exact repeat
    # - "{m,n}" or "(m,n)" -> range repeat
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

function _mask_from_char(char::Char, options::ComparisonOptions)
    # Matching is case-insensitive at parse time.
    aa = uppercase(char)
    idx = get(options.alphabet_index, aa, 0)
    idx == 0 && throw(ArgumentError("Unsupported residue '$char' for selected alphabet."))
    # The first residue in alphabet uses bit 0, second uses bit 1, etc.
    return UInt32(1) << (idx - 1)
end

function _class_mask(raw::AbstractString, options::ComparisonOptions)
    isempty(raw) && throw(ArgumentError("Empty character class is not allowed."))
    # `[^...]` syntax means "all residues except listed residues".
    invert = startswith(raw, "^")
    body = invert ? raw[nextind(raw, firstindex(raw)):end] : raw
    isempty(body) && throw(ArgumentError("Empty negated character class is not allowed."))
    mask = UInt32(0)
    for char in body
        # Class set is the union of member residue bits.
        mask |= _mask_from_char(char, options)
    end
    if invert
        # Complement inside alphabet domain only.
        mask = options.alphabet_mask & ~mask
    end
    mask == 0 && throw(ArgumentError("Character class resolves to an empty set."))
    return mask
end

function _mask_to_chars(mask::UInt32, options::ComparisonOptions; as_lowercase::Bool = false)
    chars = Char[]
    for (i, aa) in enumerate(options.alphabet)
        # Emit residues in canonical alphabet order for deterministic normalization.
        if (mask & (UInt32(1) << (i - 1))) != 0
            push!(chars, as_lowercase ? Base.lowercase(aa) : aa)
        end
    end
    return chars
end

function _mask_to_symbol(mask::UInt32, options::ComparisonOptions;
        as_lowercase::Bool = false, wildcard_symbol::String = "x")
    # Full mask is represented as wildcard, not as long explicit class.
    if mask == options.alphabet_mask
        return as_lowercase ? Base.lowercase(wildcard_symbol) : wildcard_symbol
    end
    chars = _mask_to_chars(mask, options; as_lowercase = as_lowercase)
    length(chars) == 1 && return string(chars[1])
    return "[" * join(chars) * "]"
end

function _canonical_token(position::_Position, options::ComparisonOptions)
    # Keep termini canonicalized exactly as anchors.
    if position.kind == _NTERMINUS
        return "^"
    elseif position.kind == _CTERMINUS
        return "\$"
    end
    # Residues/classes canonicalize through mask representation.
    return _mask_to_symbol(position.mask, options; as_lowercase = false, wildcard_symbol = "x")
end

# -----------------------------------------------------------------------------
# Motif parser
# -----------------------------------------------------------------------------
# The parser keeps syntax handling deterministic:
# - fixed residues, wildcard tokens (`x`, `X`, `.`), classes, negated classes,
#   and termini are parsed into `_Position` tokens;
# - optional quantifiers `{n}`, `{m,n}`, `(n)`, `(m,n)` are attached to the
#   preceding token;
# - canonical text is rebuilt from masks and quantifiers to normalize input.

function _parse_motif(motif::AbstractString, options::ComparisonOptions)
    stripped = strip(motif)
    isempty(stripped) && throw(ArgumentError("Motif cannot be empty."))

    tokens = _Token[]
    i = firstindex(stripped)
    while i <= lastindex(stripped)
        char = stripped[i]
        if isspace(char)
            # Ignore whitespace so users can provide readable motifs.
            i = nextind(stripped, i)
            continue
        end

        # Parse a single symbolic unit (position/anchor/class/wildcard).
        position = if char == '^'
            _Position(_NTERMINUS, 0)
        elseif char == '$'
            _Position(_CTERMINUS, 0)
        elseif char == 'x' || char == 'X' || char == '.'
            # Wildcards match any residue in the selected alphabet.
            _Position(_RESIDUE, options.alphabet_mask)
        elseif char == '['
            # Class token includes everything until the first closing bracket.
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
            # Anchors are single logical positions, never repeatable.
            throw(ArgumentError("Repeat quantifiers are not valid for termini in motif: $motif"))
        end
        if repeat_max == 0
            # Quantifier like {0} removes this token from every expanded variant.
            i = next_i
            continue
        end

        canonical = _canonical_token(position, options)
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
    normalized = join(getfield.(tokens, :canonical))
    return _ParsedMotif(String(motif), normalized, tokens)
end

function _variant_count(tokens::Vector{_Token})
    # Product of each token's repeat choice count, computed in BigInt to avoid overflow.
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

    # Depth-first expansion over repeat ranges.
    variants = _MotifVariant[]
    positions = _Position[]
    symbols = String[]

    function rec(ti::Int)
        if ti > length(parsed.tokens)
            isempty(positions) && return
            # Compute total motif information from concrete expanded positions.
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
            # Append repeated token positions for this branch, recurse, then backtrack.
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
