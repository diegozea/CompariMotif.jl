```@meta
CurrentModule = Main.CompariMotif
DocTestSetup = quote
    using CompariMotif
    using DataFrames
end
```

# External API

The public API is intentionally small: configure a comparison with 
[`ComparisonOptions`](@ref), compare motifs with [`compare`](@ref), and turn the results 
into a table with [`to_column_table`](@ref). These functions are built on a clean-room, 
black-box reproduction of the CompariMotif algorithm described by 
Edwards, Davey, and Shields (2008), _"CompariMotif: quick and easy comparisons 
of sequence motifs"_, [Bioinformatics](https://doi.org/10.1093/bioinformatics/btn105).

CompariMotif compares one motif against another and scores how well their
positions overlap. The best overlap is returned as a [`ComparisonResult`](@ref). If no 
significant overlap is found, the `matched` field of the result is `false`.

## Quick Start

Start by loading the package:

```julia
using CompariMotif
```

Then try a single pairwise comparison to see the main workflow. Create
[`ComparisonOptions`](@ref), call [`compare`](@ref), and inspect the returned
[`ComparisonResult`](@ref).

```jldoctest
julia> options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
ComparisonOptions(
  alphabet                = CompariMotif.ProteinAlphabet()
  min_shared_positions    = 1
  normalized_ic_cutoff    = 0.0
  matchfix                = CompariMotif.MatchFixNone
  mismatches              = 0
  allow_ambiguous_overlap = true
  max_variants            = 10000
)

julia> result = compare("RKLI", "R[KR]L[IV]", options)
ComparisonResult(
  query               = RKLI
  search              = R[KR]L[IV]
  normalized_query    = RKLI
  normalized_search   = R[RK]L[IV]
  matched             = true
  query_relationship  = Variant Match
  search_relationship = Degenerate Match
  matched_pattern     = R[rk]L[iv]
  matched_positions   = 4
  match_ic            = 3.537243573680481
  normalized_ic       = 1.0
  core_ic             = 0.8843108934201203
  score               = 4.0
  query_information   = 4.0
  search_information  = 3.537243573680481
)

julia> result.matched
true

julia> result.query_relationship
"Variant Match"
```
You can find more advanced examples of how to use this package in the 
[`FAQ / How-To`](@ref) section of the documentation.

## Matrix Comparisons

When you have many motifs, you can compare the entire collection in a single call. The 
result is a square matrix of all-vs-all pairwise comparisons, which can be converted into 
a column-oriented table for storage or downstream analysis. The table can be easily 
converted into a `DataFrame` from the `DataFrames` package or saved into a comma-separated 
values (CSV) file with the `CSV` package.

```jldoctest
julia> motifs = ["RKLI", "R[KR]L[IV]", "RxLE"];

julia> options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> results = compare(motifs, options);

julia> size(results)
(3, 3)

julia> table = to_column_table(results);

julia> length(table.query)
9

```

## Canonicalization

[`normalize_motif`](@ref) parses a motif and renders it in the deterministic syntax 
used internally by the comparison pipeline. You do not need to call this function directly.

```jldoctest
julia> normalize_motif("r[kR].{0,1}l")
"R[RK]x{0,1}L"
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
