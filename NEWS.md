## CompariMotif.jl Release Notes

### Changes from 0.1.0 to main

- *[Breaking change]* Alphabet selectors now use marker objects instead of
  symbols or strings in the public API.
  - Use `ProteinAlphabet()`, `DNAAlphabet()`, or `RNAAlphabet()` with
  `ComparisonOptions` and `normalize_motif`.
  - Migration examples:
    - `ComparisonOptions(; alphabet = :dna)` becomes
      `ComparisonOptions(; alphabet = DNAAlphabet())`
    - `normalize_motif("AUG"; alphabet = :rna)` becomes
      `normalize_motif("AUG"; alphabet = RNAAlphabet())`
- *[Breaking change]* `ComparisonOptions` now keeps the selected alphabet marker in the public
  struct and display output instead of exposing derived internal alphabet
  caches.
- *[Breaking change]* `ComparisonOptions.matchfix` now uses symbols instead of
  the exported `MatchFixMode` enum and enum values.
  - Valid values are exactly `:none`, `:query_fixed`, `:search_fixed`, and `:both_fixed`.
  - Strings and legacy enum names are no longer accepted.
- `compare` now also accepts a single query motif as the first argument and 
  a  collection of target motifs as the second argument.
- Added support for non-uniform residue frequencies through
  `ComparisonOptions(; residue_frequencies = ...)`, so information-content
  scoring can use a custom DNA, RNA, or protein background model instead of a
  uniform distribution. The new option validates alphabet coverage and is
  covered by dedicated oracle-backed regression tests.
- Fixed the comparison pipeline to perform the paper's precise exact-match and
  exact-subsequence pre-pass before the full sliding-window overlap search.
  This improves candidate selection for motifs with alternation and repeat
  expansion.
- Added a multi-page documentation manual with separate sections for the
  external API, the internal API and comparison pipeline, and a FAQ/how-to
  guide.
- Expanded paper-backed regression coverage to Figures 1, 2, and 3, and added
  a tested internal pipeline walkthrough that stays aligned with the manual.
