# Canonical alphabets used by the parser and IC calculations.
# The order is stable and intentionally used as bit-index order in residue masks.
const _PROTEIN_ALPHABET = collect("ARNDCQEGHILKMFPSTWYV")
const _DNA_ALPHABET = collect("ACGT")
const _RNA_ALPHABET = collect("ACGU")

# Shared docstring snippets to keep user-facing documentation consistent
# across `normalize_motif`, `compare`, and `to_column_table`.
const _DOC_OPTIONS_REF = "Configure thresholds and matching semantics with [`ComparisonOptions`](@ref)."
const _DOC_RESULT_REF = "Returns a [`ComparisonResult`](@ref)."
const _DOC_NORMALIZE_REF = "Use [`normalize_motif`](@ref) for deterministic motif canonicalization."
const _DOC_COMPARE_REF = "Compute similarities with [`compare`](@ref)."
const _DOC_TABLE_REF = "Convert results to column tables with [`to_column_table`](@ref)."
const _DOC_VARIANT_SIZE_REF = "The result matrix has size `(length(motifs), length(db))`."

# Relationship render tables (human-readable words).
const _RELATIONSHIP_TYPE_WORDS = ("Exact", "Variant", "Degenerate", "Complex")
const _RELATIONSHIP_LENGTH_WORDS = ("Match", "Parent", "Subsequence", "Overlap")
