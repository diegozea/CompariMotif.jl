## CompariMotif.jl Release Notes

### Changes from 0.1.0 to 0.2.0

- *[Breaking change]* Alphabet selectors now use marker objects instead of
  symbols or strings in the public API.
  - Use `ProteinAlphabet()`, `DNAAlphabet()`, or `RNAAlphabet()` with
  `ComparisonOptions`.
  - Migration examples:
    - `ComparisonOptions(; alphabet = :dna)` becomes
      `ComparisonOptions(; alphabet = DNAAlphabet())`
- *[Breaking change]* `normalize_motif` is no longer part of the public API.
  - Canonicalized motif strings remain available through
    `ComparisonResult.normalized_query`, `ComparisonResult.normalized_search`,
    and `to_column_table`.
- *[Breaking change]* `ComparisonOptions` now keeps the selected alphabet marker in the public
  struct and display output instead of exposing derived internal alphabet
  caches.
- *[Breaking change]* `ComparisonOptions.matchfix` now uses symbols instead of
  the exported `MatchFixMode` enum and enum values.
  - Valid values are exactly `:none`, `:query_fixed`, `:search_fixed`, and `:both_fixed`.
  - Strings and legacy enum names are no longer accepted.
- `compare` now also accepts a single query motif as the first argument and 
  a  collection of target motifs as the second argument.
- When multiple candidate overlaps tie on `match_ic` and matched positions,
  `compare` now prefers the higher final `score`, and otherwise keeps the
  first candidate encountered in the oracle-backed shift/branch scan order.
  This improves deterministic branch selection while matching the current
  tie-break fixtures.
- Added support for non-uniform residue frequencies through
  `ComparisonOptions(; residue_frequencies = ...)`, so information-content
  scoring can use a custom DNA, RNA, or protein background model instead of a
  uniform distribution. The new option validates alphabet coverage and is
  covered by dedicated oracle-backed regression tests.
- Fixed the comparison pipeline to perform the paper's precise exact-match and
  exact-subsequence pre-pass before the full sliding-window overlap search.
  This improves candidate selection for motifs with alternation and repeat
  expansion.
- Fixed the final shift-selection tie-break to follow the shift order inferred
  from black-box oracle tie cases when `match_ic`, matched positions, and
  `score` are all tied.
- Fixed alternation-heavy motifs to enforce `max_variants` during grouped-branch
  parsing, so pathological branch growth is rejected before full parse-time
  expansion.
- Wildcard aliases `x`, `X`, and `.` are now treated identically in canonical
  parsing and normalization, even in top-level grouped alternations. This is an
  intentional divergence from the current oracle quirk, where top-level
  `(Q|.)`-style grouped branches can disappear entirely while `(Q|x)` and
  `(Q|X)` retain the non-wildcard branch and only lose the wildcard-alias
  branch.
- Positive character classes are intentionally treated as sets, so duplicate
  residues like `[AA]` and `[A]` compare identically even though the current
  oracle may score them differently.
- Added a multi-page documentation manual with separate sections for the
  external API, the internal API and comparison pipeline, and a FAQ/how-to
  guide, plus a dedicated regex syntax reference for supported parser edge
  cases.
- Expanded paper-backed regression coverage to Figures 1, 2, and 3, and added
  a tested internal pipeline walkthrough that stays aligned with the manual.
