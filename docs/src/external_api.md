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
MatchFixMode
ComparisonOptions
ComparisonResult
normalize_motif
compare
to_column_table
```
