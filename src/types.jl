"""
    ComparisonResult

Result record produced by [`compare`](@ref) for one query/search motif pair.

Fields:
- `query`, `search`: original input motifs.
- `normalized_query`, `normalized_search`: canonicalized motifs used internally.
- `matched`: whether the best-scoring valid alignment passed all thresholds.
- `query_relationship`, `search_relationship`: directional two-word labels.
  The first word describes specificity (`Exact`, `Variant`, `Degenerate`,
  `Complex`); the second describes coverage (`Match`, `Parent`,
  `Subsequence`, `Overlap`).
- `matched_pattern`: compact rendering of the selected overlap; lowercase
  symbols mark positions broadened by ambiguity or wildcard handling.
- `matched_positions`: number of informative aligned positions; dual wildcards
  are excluded.
- `match_ic`: raw information content captured by the selected alignment.
- `normalized_ic`: `match_ic` scaled by the less informative of the two motif
  variants, making hits easier to compare across motif lengths and specificity.
- `core_ic`: fraction of aligned-core information retained by the selected
  alignment (`match_ic / sum(max(position_ic(query_i), position_ic(search_i)))`
  across non-dual-wildcard aligned positions).
- `score`: derived summary score (`normalized_ic * matched_positions`).
- `query_information`, `search_information`: total information content for the
  winning query and search motif variants.

The relationship fields are asymmetric by design, so one side of the same hit
can read `Variant Subsequence` while the other reads `Degenerate Parent`.
When `matched == false`, the relationship labels are `No Match`,
`matched_pattern` is empty, and all position/score/information totals stay at
their zero defaults.

See also [`ComparisonOptions`](@ref), [`compare`](@ref), [`to_column_table`](@ref).
"""
Base.@kwdef struct ComparisonResult
    query::String
    search::String
    normalized_query::String
    normalized_search::String
    matched::Bool = false
    query_relationship::String = "No Match"
    search_relationship::String = "No Match"
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

# Internal position kind:
# - residue: position with a residue set mask,
# - termini: boundary markers (`^` or `$`) that are matched as anchors.
@enum _PositionKind::UInt8 begin
    _RESIDUE = 0
    _NTERMINUS = 1
    _CTERMINUS = 2
end

"""
    ResidueClass

Residue set encoded as a [`ResidueMask`](@ref).
"""
struct ResidueClass
    mask::ResidueMask
end

# Core atom used everywhere in comparison code.
# `mask` is meaningful only when `kind == _RESIDUE`.
struct _Position
    kind::_PositionKind
    mask::ResidueMask
end

# Parser token = one parsed position plus optional repeat range.
# `canonical` stores the deterministic textual representation for normalization.
struct _Token
    position::_Position
    min_repeat::Int
    max_repeat::Int
    canonical::String
end

# Parsed motif keeps:
# - original user input,
# - normalized canonical text,
# - token sequence with possible repeat ranges,
# - branch alternatives for `(motif1)|(motif2)`-style motifs.
"""
    _ParsedMotif

Internal parsed representation of one user-supplied motif.

Fields:
- `original`: motif text exactly as supplied by the caller.
- `normalized`: canonical motif text used for deterministic comparisons.
- `tokens`: token sequence for the first parsed branch.
- `alternatives`: token sequence for every expanded top-level alternation branch.
"""
struct _ParsedMotif
    original::String
    normalized::String
    tokens::Vector{_Token}
    alternatives::Vector{Vector{_Token}}
end

# Expanded motif variant after resolving repeat ranges.
# `positions` is the concrete position sequence used for alignment.
"""
    _MotifVariant

Concrete motif variant obtained after expanding bounded repeat ranges.

Fields:
- `positions`: fixed sequence of parsed positions used during alignment.
- `normalized`: canonical motif text for this expanded variant.
- `information`: total information content of the variant.
"""
struct _MotifVariant
    positions::Vector{_Position}
    normalized::String
    information::Float64
end

"""
    ComparisonOptions

Reusable configuration object for CompariMotif comparisons.

Construct once with [`ComparisonOptions(; kwargs...)`](@ref) and reuse across
many [`compare`](@ref) calls.

# Keywords
- `alphabet = ProteinAlphabet()`: comparison alphabet (`ProteinAlphabet()`,
  `DNAAlphabet()`, or `RNAAlphabet()`).
- `residue_frequencies::Union{Nothing, AbstractDict{Char,<:Real}} = nothing`:
  optional background residue frequencies for information-content scoring.
  When omitted, CompariMotif uses a uniform frequency distribution. Provided
  dictionaries must define every residue in the selected alphabet, use
  strictly positive finite values. Frequencies are normalized internally to sum to `1.0`.
- `min_shared_positions::Int = 2`: minimum number of matched, non-wildcard
  positions required for a hit.
- `normalized_ic_cutoff::Real = 0.5`: minimum normalized information content.
- `matchfix::Symbol = :none`: fixed-position matching mode. Accepted values are
  exactly `:none`, `:query_fixed`, `:search_fixed`, and `:both_fixed`.
- `mismatches::Int = 0`: tolerated count of defined-position mismatches.
- `allow_ambiguous_overlap::Bool = true`: whether partial class overlaps are
  allowed as complex matches.
- `max_variants::Int = 10_000`: maximum expanded variants per motif.

# Examples
```jldoctest
julia> using CompariMotif

julia> options = ComparisonOptions();

julia> options.alphabet isa ProteinAlphabet
true

julia> options.matchfix == :none
true
```

See also [`compare`](@ref), [`ComparisonResult`](@ref).
"""
struct ComparisonOptions
    # Public alphabet selector used to reconstruct or display the configuration.
    alphabet::_AlphabetValue
    # Optional normalized background frequencies keyed by canonical alphabet residue.
    residue_frequencies::Union{Nothing, Dict{Char, Float64}}
    # Minimum matched non-wildcard positions for a valid hit.
    min_shared_positions::Int
    # Minimum normalized information score for a valid hit.
    normalized_ic_cutoff::Float64
    # Fixed-position matching policy: :none, :query_fixed, :search_fixed or :both_fixed
    matchfix::Symbol
    # Allowed number of mismatch positions inside an overlap.
    mismatches::Int
    # Whether partial overlap of two ambiguous sets is allowed.
    allow_ambiguous_overlap::Bool
    # Safety limit for repeat expansion combinatorics.
    max_variants::Int

    function ComparisonOptions(
            alphabet::_AlphabetValue,
            residue_frequencies::Union{Nothing, Dict{Char, Float64}},
            min_shared_positions::Int,
            normalized_ic_cutoff::Float64,
            matchfix::Symbol,
            mismatches::Int,
            allow_ambiguous_overlap::Bool,
            max_variants::Int
    )
        return new(
            alphabet,
            residue_frequencies,
            min_shared_positions,
            normalized_ic_cutoff,
            _coerce_matchfix(matchfix),
            mismatches,
            allow_ambiguous_overlap,
            max_variants
        )
    end
end

@def_pprint mime_types="text/plain" base_show=true ComparisonOptions

@inline _alphabet_spec(options::ComparisonOptions) = _alphabet_spec(options.alphabet)

"""
    overlaps(a::ResidueClass, b::ResidueClass)::Bool

Return `true` when two residue classes share at least one residue.
"""
@inline overlaps(a::ResidueClass, b::ResidueClass) = !iszero(a.mask & b.mask)

"""
    unionclass(a::ResidueClass, b::ResidueClass)::ResidueClass

Return the set-union of two residue classes.
"""
@inline unionclass(a::ResidueClass, b::ResidueClass) = ResidueClass(a.mask | b.mask)

"""
    is_subset(a::ResidueClass, b::ResidueClass)::Bool

Return `true` when every residue in `a` is also in `b`.
"""
@inline is_subset(a::ResidueClass, b::ResidueClass) = iszero(a.mask & ~b.mask)

"""
    is_wildcard(a::ResidueClass, opts::ComparisonOptions)::Bool

Return `true` when the residue class spans the full selected alphabet.
"""
@inline is_wildcard(a::ResidueClass, opts::ComparisonOptions) = a.mask ==
                                                                _alphabet_spec(opts.alphabet).mask

"""
    is_fixed(a::ResidueClass)::Bool

Return `true` when the residue class contains exactly one residue.
"""
@inline is_fixed(a::ResidueClass) = count_ones(a.mask) == 1

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

# A scored candidate alignment between one query variant and one search variant.
# This carries enough detail to build the final public `ComparisonResult`.
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
