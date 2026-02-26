# Canonical alphabets used by the parser and IC calculations.
# The order is stable and intentionally used as bit-index order in residue masks.
const _PROTEIN_ALPHABET = collect("ACDEFGHIKLMNPQRSTVWY")
const _DNA_ALPHABET = collect("ACGT")

# Shared docstring snippets to keep user-facing documentation consistent
# across `normalize_motif`, `compare`, and `write_results_tsv`.
const _DOC_OPTIONS_REF = "Configure thresholds and matching semantics with [`ComparisonOptions`](@ref)."
const _DOC_RESULT_REF = "Returns a [`ComparisonResult`](@ref)."
const _DOC_NORMALIZE_REF = "Use [`normalize_motif`](@ref) for deterministic motif canonicalization."
const _DOC_COMPARE_REF = "Compute similarities with [`compare`](@ref)."
const _DOC_TSV_REF = "Persist matrices with [`write_results_tsv`](@ref)."
const _DOC_VARIANT_SIZE_REF = "The result matrix has size `(length(motifs), length(db))`."

# Relationship render tables:
# - words for human-readable result fields,
# - compact codes for regression fixtures and TSV output.
const _RELATIONSHIP_TYPE_WORDS = ("Exact", "Variant", "Degenerate", "Complex")
const _RELATIONSHIP_TYPE_CODES = ("e", "v", "d", "c")
const _RELATIONSHIP_LENGTH_WORDS = ("Match", "Parent", "Subsequence", "Overlap")
const _RELATIONSHIP_LENGTH_CODES = ("m", "p", "s", "o")
