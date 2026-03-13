```@meta
CurrentModule = Main.CompariMotif
```

# External API

The public API is intentionally small: configure a comparison with 
[`ComparisonOptions`](@ref), compare motifs with [`compare`](@ref), and turn the results 
into a table with [`to_column_table`](@ref). These functions are built on a clean-room, 
black-box reproduction of the CompariMotif algorithm described by
[Edwards2008CompariMotif](@citet).

CompariMotif compares one motif against another and scores how well their
positions overlap. The best overlap is returned as a [`ComparisonResult`](@ref). If no 
significant overlap is found, the `matched` field of the result is `false`.

## Quick Start

Start by loading the package:

```@repl external_api_examples
using CompariMotif
```

Then try a single pairwise comparison to see the main workflow. Create
[`ComparisonOptions`](@ref), call [`compare`](@ref), and inspect the returned
[`ComparisonResult`](@ref).

```@repl external_api_examples
options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
result = compare("RKLI", "R[KR]L[IV]", options)
```

## Interpreting `ComparisonResult`

Each [`ComparisonResult`](@ref) summarizes the best accepted overlap between one
query motif and one search motif. The relationship fields are directional:
`query_relationship` uses the query as the reference point, while
`search_relationship` describes the same alignment from the search side. That
is why `Variant` and `Degenerate`, and `Parent` and `Subsequence`, normally
appear as complementary pairs. For example, one side of a hit can read
`Degenerate Parent` while the other reads `Variant Subsequence`.

The first word in the relationship label explains how residue choices compare
along the overlap: `Exact` means the aligned residue sets coincide, `Variant`
means the query is narrower, `Degenerate` means the query is broader, and
`Complex` means the overlap mixes those cases or uses partly overlapping
ambiguous classes. The second word explains coverage: `Match` for full-length
coverage on both motifs, `Parent` when the query contains the search,
`Subsequence` when the query is contained in the search, and `Overlap` when the
best hit uses only part of each motif. The nomenclature used here is identical to 
that of [Edwards2008CompariMotif](@citet), which is nicely summarized in _Figure 2_.

The score fields then tell you how strong that overlap is. `match_ic` is the
raw information captured by the best alignment, `normalized_ic` puts that value
on a comparable scale across motifs of different specificity, `core_ic` is the
average information per aligned core position, and `score` combines
`normalized_ic` with `matched_positions` for ranking. `matched_positions`
counts informative aligned positions only, so positions where both motifs use a
wildcard do not contribute to that total. `matched_pattern` is a compact
rendering of the winning overlap: uppercase symbols usually mark clean exact
agreement, whereas lowercase symbols mark positions broadened by ambiguity or
wildcard handling.

The same fields become table columns when you call [`to_column_table`](@ref).
In that tabular form, `query` and `search` keep the original motif strings,
while `normalized_query` and `normalized_search` expose the canonical forms
used internally during comparison. Unmatched results keep `matched = false`,
use `No Match` relationship labels, leave `matched_pattern` empty, and set
`matched_positions`, all score fields, and both information totals to `0` or
`0.0`.

You can find more advanced examples of how to use this package in the
[FAQ / How-To](@ref) section of the documentation.

## Non-Uniform Residue Frequencies

By default, information content uses a uniform residue frequency distribution.
To score motifs against a custom background model, pass
`residue_frequencies = Dict{Char,Float64}(...)` when constructing
[`ComparisonOptions`](@ref).

```@repl external_api_examples
dna_freqs = Dict('A' => 0.3, 'C' => 0.2, 'G' => 0.2, 'T' => 0.3)
weighted = ComparisonOptions(;
    alphabet = DNAAlphabet(),
    residue_frequencies = dna_freqs,
    min_shared_positions = 1,
    normalized_ic_cutoff = 0.0,
)
compare("ATG", "[AGT]TG", weighted)
```

## Matrix Comparisons

When you have many motifs, you can compare the entire collection in a single call. The 
result is a square matrix of all-vs-all pairwise comparisons, which can be converted into 
a column-oriented table for storage or downstream analysis. The table can be easily 
converted into a `DataFrame` from the `DataFrames` package or saved into a comma-separated 
values (CSV) file with the `CSV` package.

```@repl external_api_examples
motifs = ["RKLI", "R[KR]L[IV]", "RxLE"]
results = compare(motifs, options);
size(results)
table = to_column_table(results)
```

When you want to search a single query motif against a database of targets,
pass the query as a string and the targets as a vector. This is useful when you
already know the motif you want to look up and want to inspect or export all
hits against a target set.

```@repl external_api_examples
query_hits = compare("RKLI", ["RKLI", "RxLE"], options);
query_hits[1, 1].matched
query_hits[1, 2].matched
```

## Canonicalization

[`normalize_motif`](@ref) parses a motif and renders it in the deterministic syntax 
used internally by the comparison pipeline. You do not need to call this function directly.

```@repl external_api_examples
normalize_motif("r[kR].{0,1}l")
```

## Reference

The full public API is listed below. Use these docstrings when you need the
precise meaning of an option, result field, or helper function.

```@docs
ProteinAlphabet
DNAAlphabet
RNAAlphabet
ComparisonOptions
ComparisonResult
normalize_motif
compare
to_column_table
```
