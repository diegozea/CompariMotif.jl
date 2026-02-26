"""
    ComparisonResult

Result record produced by [`compare`](@ref) for one query/search motif pair.

Fields:
- `query`, `search`: original input motifs.
- `normalized_query`, `normalized_search`: canonicalized motifs used internally.
- `matched`: whether the best-scoring valid alignment passed all thresholds.
- `query_relationship`, `search_relationship`: human-readable relationship labels.
- `query_relationship_code`, `search_relationship_code`: compact relationship codes.
- `matched_pattern`: consensus/overlap pattern for the selected alignment.
- `matched_positions`: count of matched non-wildcard positions.
- `match_ic`: total information content for matched positions.
- `normalized_ic`: `match_ic` normalized by the lower motif information content.
- `core_ic`: information content normalized by core overlap length.
- `score`: derived summary score (`normalized_ic * matched_positions`).
- `query_information`, `search_information`: total information content per motif.

See also [`ComparisonOptions`](@ref), [`normalize_motif`](@ref), [`write_results_tsv`](@ref).
"""
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

# Internal position kind:
# - residue: position with a residue set mask,
# - termini: boundary markers (`^` or `$`) that are matched as anchors.
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

Used by the `matchfix` keyword in [`ComparisonOptions`](@ref).
"""
@enum MatchFixMode::UInt8 begin
    MatchFixNone = 0
    MatchFixQueryFixed = 1
    MatchFixSearchFixed = 2
    MatchFixBothFixed = 3
end

# Core atom used everywhere in comparison code.
# `mask` is meaningful only when `kind == _RESIDUE`.
struct _Position
    kind::_PositionKind
    mask::UInt32
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
# - token sequence with possible repeat ranges.
struct _ParsedMotif
    original::String
    normalized::String
    tokens::Vector{_Token}
end

# Expanded motif variant after resolving repeat ranges.
# `positions` is the concrete position sequence used for alignment.
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

See also [`MatchFixMode`](@ref), [`compare`](@ref), [`ComparisonResult`](@ref).
"""
struct ComparisonOptions
    # Ordered alphabet used for mask generation and IC scaling.
    alphabet::Vector{Char}
    # Maps residue character -> 1-based bit position.
    alphabet_index::Dict{Char, Int}
    # Bitmask with all alphabet residues enabled (wildcard mask).
    alphabet_mask::UInt32
    # Precomputed log base for IC normalization (`log(N)`).
    log_base::Float64
    # Minimum matched non-wildcard positions for a valid hit.
    min_shared_positions::Int
    # Minimum normalized information score for a valid hit.
    normalized_ic_cutoff::Float64
    # Fixed-position matching policy.
    matchfix::MatchFixMode
    # Allowed number of mismatch positions inside an overlap.
    mismatches::Int
    # Whether partial overlap of two ambiguous sets is allowed.
    allow_ambiguous_overlap::Bool
    # Safety limit for repeat expansion combinatorics.
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
