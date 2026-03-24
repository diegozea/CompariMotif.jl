const _MATCHFIX_ERROR_MESSAGE = "`matchfix` must be one of :none, :query_fixed, :search_fixed, :both_fixed."

_matchfix_argument_error() = ArgumentError(_MATCHFIX_ERROR_MESSAGE)

"""
    _coerce_matchfix(mode::Symbol)::Symbol

Validate `matchfix` and return one of the canonical symbols accepted by
[`ComparisonOptions`](@ref).
"""
function _coerce_matchfix(mode::Symbol)
    # Keep only the public canonical symbols so storage and comparisons stay simple.
    if mode === :none
        return :none
    elseif mode === :query_fixed
        return :query_fixed
    elseif mode === :search_fixed
        return :search_fixed
    elseif mode === :both_fixed
        return :both_fixed
    end
    throw(_matchfix_argument_error())
end

function _canonicalize_residue_frequencies(::Nothing, ::_AlphabetSpec)
    return nothing
end

"""
    _canonicalize_residue_frequencies(residue_frequencies, spec::_AlphabetSpec)

Validate and normalize user-supplied residue frequencies into canonical
alphabet keys whose probabilities sum to `1.0`.

# Examples
```jldoctest
julia> using CompariMotif

julia> spec = CompariMotif._alphabet_spec(CompariMotif.DNAAlphabet());

julia> freqs = CompariMotif._canonicalize_residue_frequencies(
           Dict('a' => 7.0, 'c' => 2.0, 'g' => 1.0, 't' => 1.0),
           spec
       );

julia> isapprox([freqs[aa] for aa in spec.chars], [7, 2, 1, 1] / 11)
true
```
"""
function _canonicalize_residue_frequencies(
        residue_frequencies::AbstractDict{Char, <:Real},
        spec::_AlphabetSpec
)
    normalized = Dict{Char, Float64}()
    for (key, value) in pairs(residue_frequencies)
        aa = uppercase(key)
        haskey(spec.index, aa) || throw(ArgumentError(
            "`residue_frequencies` contains unsupported residue '$key' for the selected alphabet."
        ))
        haskey(normalized, aa) && throw(ArgumentError(
            "`residue_frequencies` contains duplicate residue '$aa' after case normalization."
        ))
        numeric = Float64(value)
        isfinite(numeric) || throw(
            ArgumentError("`residue_frequencies` values must be finite.")
        )
        numeric > 0.0 || throw(
            ArgumentError(
            "`residue_frequencies` values must be positive; use pseudocounts."
        ))
        normalized[aa] = numeric
    end

    missing_residues = Char[]
    for aa in spec.chars
        haskey(normalized, aa) || push!(missing_residues, aa)
    end
    isempty(missing_residues) || throw(ArgumentError(
        "`residue_frequencies` must define every residue in the selected alphabet; missing: " *
        join(missing_residues, ", ")
    ))

    # Scale before summation so extremely large or uneven user weights normalize
    # without overflowing the Float64 accumulator.
    scale = maximum(values(normalized))
    scale > 0.0 || throw(
        ArgumentError("`residue_frequencies` must have positive total mass.")
    )
    total = 0.0
    for value in values(normalized)
        total += value / scale
    end
    isfinite(total) && total > 0.0 || throw(
        ArgumentError("`residue_frequencies` must normalize to a finite positive total mass.")
    )

    # Deterministic downstream access uses `spec.chars`, not dictionary
    # iteration order, so a plain `Dict` is sufficient here.
    canonical_frequencies = Dict{Char, Float64}()
    for aa in spec.chars
        probability = (normalized[aa] / scale) / total
        probability > 0.0 || throw(ArgumentError(
            "`residue_frequencies` span too large a dynamic range to normalize safely in Float64."
        ))
        canonical_frequencies[aa] = probability
    end
    return canonical_frequencies
end

function _uniform_residue_frequency_vector(spec::_AlphabetSpec)
    n = length(spec.chars)
    return fill(inv(n), n)
end

function _residue_frequency_vector(options::ComparisonOptions, spec::_AlphabetSpec)
    frequencies = options.residue_frequencies
    frequencies === nothing && return _uniform_residue_frequency_vector(spec)
    return [frequencies[aa] for aa in spec.chars]
end

# Keep positional and keyword construction on the same validation/canonicalization path.
function _build_comparison_options(
        alphabet::_AlphabetValue,
        residue_frequencies::Union{Nothing, AbstractDict{Char, <:Real}},
        min_shared_positions::Int,
        normalized_ic_cutoff::Real,
        matchfix::Symbol,
        mismatches::Int,
        allow_ambiguous_overlap::Bool,
        max_variants::Int
)
    spec = _alphabet_spec(alphabet)
    normalized_frequencies = _canonicalize_residue_frequencies(residue_frequencies, spec)
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

    return ComparisonOptions(
        alphabet,
        normalized_frequencies,
        min_shared_positions,
        Float64(normalized_ic_cutoff),
        matchfix,
        mismatches,
        allow_ambiguous_overlap,
        max_variants
    )
end

"""
    ComparisonOptions(; kwargs...)::ComparisonOptions

Construct a reusable options object for motif comparisons.

```jldoctest
julia> using CompariMotif

julia> ComparisonOptions(; alphabet = DNAAlphabet()).alphabet isa DNAAlphabet
true
```
"""
function ComparisonOptions(;
        alphabet::_AlphabetValue = ProteinAlphabet(),
        residue_frequencies::Union{Nothing, AbstractDict{Char, <:Real}} = nothing,
        min_shared_positions::Int = 2,
        normalized_ic_cutoff::Real = 0.5,
        matchfix::Symbol = :none,
        mismatches::Int = 0,
        allow_ambiguous_overlap::Bool = true,
        max_variants::Int = 10_000
)
    return _build_comparison_options(
        alphabet,
        residue_frequencies,
        min_shared_positions,
        normalized_ic_cutoff,
        matchfix,
        mismatches,
        allow_ambiguous_overlap,
        max_variants
    )
end

function ComparisonOptions(
        alphabet::_AlphabetValue,
        residue_frequencies::Union{Nothing, AbstractDict{Char, <:Real}},
        min_shared_positions::Int,
        normalized_ic_cutoff::Real,
        matchfix::Symbol,
        mismatches::Int,
        allow_ambiguous_overlap::Bool,
        max_variants::Int
)
    return _build_comparison_options(
        alphabet,
        residue_frequencies,
        min_shared_positions,
        normalized_ic_cutoff,
        matchfix,
        mismatches,
        allow_ambiguous_overlap,
        max_variants
    )
end
