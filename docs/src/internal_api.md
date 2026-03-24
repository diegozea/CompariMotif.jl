```@meta
CurrentModule = Main.CompariMotif
```

# Internal API & Pipeline

This page documents the current implementation pipeline. Everything here is
private and may change between releases; the stable package contract is the
[External API](external_api.md).

The running example below mirrors the worked comparison in Figure 1 of
[Edwards2008CompariMotif](@citet).

## Pipeline Overview

The current code path for a pairwise comparison is:

1. [`_parse_motif`](@ref) strips whitespace, expands grouping and alternation,
   parses tokens, and produces a canonical normalized motif string.
2. [`_expand_variants`](@ref) resolves bounded repeat ranges into concrete motif
   variants with precomputed information content.
3. [`_find_precise_match`](@ref) checks exact full-length and exact
   subsequence relationships first to seed the best candidate before the full
   overlap search.
4. [`_evaluate_alignment`](@ref) scores each overlap candidate considered
   across the expanded variant pairs and relative shifts.
5. [`_compare_parsed`](@ref) keeps the best candidate and materializes the
   public [`ComparisonResult`](@ref).

## Figure 1 Worked Example

### 1. Parse and normalize the motifs

Parsing is the syntax-level normalization step. It canonicalizes residue
classes and wildcard notation, records any alternation branches, and preserves
bounded repeats in normalized form so later stages can expand them
deliberately. Wildcard aliases are canonicalized to `.` even where the current
oracle has a grouped-alternation quirk, because the package treats `x`, `X`,
and `.` as intentionally equivalent syntax. Positive character classes are also
treated as sets, so duplicate residues are discarded even though the oracle can
score them differently.

```@repl internal_api_pipeline
using CompariMotif # hide
options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);
parsed_query = CompariMotif._parse_motif("[KR].L.{0,1}[FYLIVMP]", options);
parsed_search = CompariMotif._parse_motif("R.LE", options);
parsed_query.normalized
parsed_search.normalized
```

### 2. Expand concrete variants

Variant expansion converts each parsed branch into one or more concrete motif
variants with explicit positions and precomputed information content. This is
the stage where repeat ranges become enumerated sequences, so all downstream
alignment and scoring logic works with concrete variant objects rather than
quantified syntax.

```@repl internal_api_pipeline
spec = CompariMotif._alphabet_spec(options.alphabet);
query_variants = CompariMotif._expand_variants(parsed_query, options, spec);
search_variants = CompariMotif._expand_variants(parsed_search, options, spec);
[variant.normalized for variant in query_variants]
round.([variant.information for variant in query_variants], digits = 3)
only(search_variants).normalized
round(only(search_variants).information, digits = 3)
```

### 3. Check precise matches before overlap scoring

The precise-match pass looks for exact full-length and exact subsequence
relationships among the expanded variants before the broader overlap search
runs. The current implementation does not add a separate check for whether two
motifs contain enough shared amino acids in any position to merit further
comparison; after the exact-match pass, it proceeds directly to the
sliding-window overlap search. Any exact hit seeds the current best candidate,
but it still does not short-circuit later evaluation of other overlaps. We have tried 
implementing that, as it was suggested in the paper, but it did not improve 
performance in practice.

```@repl internal_api_pipeline
found_precise, best_precise = CompariMotif._find_precise_match(query_variants, search_variants, options, spec);
found_precise
```

### 4. Score the best overlap

Alignment scoring evaluates one query variant against one search variant at a
specific relative shift. Each candidate carries the matched pattern, matched
positions, relationship labels, and the information-content-derived metrics
used for ranking. The current implementation orders candidates by higher
`match_ic`, then `matched_positions`, then `score`. If all of those still tie,
the first candidate encountered in the shift scan inferred from black-box
oracle tie cases is kept.
`normalized_ic`, `core_ic`, and `score` are still materialized on the
candidate for inspection and output.

```@repl internal_api_pipeline
query_variant = query_variants[2];
search_variant = only(search_variants);
candidate = CompariMotif._evaluate_alignment(query_variant, search_variant, 0, options, spec);
candidate.matched_pattern
candidate.matched_positions
round(candidate.normalized_ic, digits = 3)
```

### 5. Materialize the public result

After all precise matches and overlap candidates have been considered,
`_compare_parsed` keeps the strongest candidate and materializes it as the
public [`ComparisonResult`](@ref). This final step copies the winning
alignment's relationships, matched pattern, and information-content summary
into the stable API object returned by [`compare`](@ref).

```@repl internal_api_pipeline
result = compare("[KR].L.{0,1}[FYLIVMP]", "R.LE", options);
(result.query_relationship, result.search_relationship)
round(result.match_ic, digits = 3)
round(result.score, digits = 3)
```

## Internal Reference

```@docs
_ParsedMotif
_MotifVariant
_parse_motif
_expand_variants
_find_precise_match
_evaluate_alignment
_compare_parsed
```
