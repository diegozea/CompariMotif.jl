# CompariMotif.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/CompariMotif.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://diegozea.github.io/CompariMotif.jl/dev)
[![CI](https://github.com/diegozea/CompariMotif.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/diegozea/CompariMotif.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/diegozea/CompariMotif.jl/badge.svg?branch=main)](https://coveralls.io/github/diegozea/CompariMotif.jl?branch=main)
[![codecov](https://codecov.io/gh/diegozea/CompariMotif.jl/graph/badge.svg?branch=main)](https://codecov.io/gh/diegozea/CompariMotif.jl?branch=main)

Clean-room, unofficial Julia implementation of the motif–motif comparison strategy 
described by Edwards, Davey and Shields (Bioinformatics 24(10):1307–1309, 2008). It 
supports the comparison of protein, DNA and RNA motifs, represented as regular 
expressions.

## API

Most workflows follow the same pattern: create a `ComparisonOptions` object
once, run one of the `compare` methods, and export results with
`to_column_table`.

### `ComparisonOptions`

- `ComparisonOptions(; kwargs...)`

Create a reusable options object that holds the alphabet, thresholds, and
matching rules for a comparison run. In practice, you normally build one
`ComparisonOptions` value at the start of an analysis and pass it to every
`compare` call in that workflow.

### `compare`

- `compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)::ComparisonResult`
- `compare(query::AbstractString, targets::AbstractVector{<:AbstractString}, options::ComparisonOptions)::Matrix{ComparisonResult}`
- `compare(motifs::AbstractVector{<:AbstractString}, db::AbstractVector{<:AbstractString}, options::ComparisonOptions)::Matrix{ComparisonResult}`
- `compare(motifs::AbstractVector{<:AbstractString}, options::ComparisonOptions)::Matrix{ComparisonResult}`

Use pairwise comparison when you want to inspect the relationship between one
known query motif and one known target motif. Use query-vs-targets mode when
you want to search a single query motif against a database of target motifs.
Use collection-vs-collection mode when you want to compare two motif sets, and
use `compare(motifs, options)` for all-vs-all comparisons within one set.

### `to_column_table`

- `to_column_table(result_or_results)::NamedTuple`

Use `to_column_table` to turn a single result or a collection of results into a
column-oriented table that can be passed directly to `DataFrame` or
`CSV.write`. This is the simplest way to persist pairwise hits,
query-vs-target search results, or larger database comparisons.

## Minimal example

```julia
using CompariMotif
using DataFrames

motifs = ["RKLI", "R[KR]L[IV]", "[KR]xLx[FYLIMVP]", "RxLE"]
options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
results = compare(motifs, options)

results[3, 4]  # single pair summary
table = to_column_table(results)
df = DataFrame(table)
```

`to_column_table` output can also be written with `CSV.write("comparimotif_results.tsv", table)`.

## Interpreting results

See [`Interpreting ComparisonResult`](https://diegozea.github.io/CompariMotif.jl/external_api/#Interpreting-ComparisonResult)
in the external API documentation for the full explanation of relationship
labels, scoring fields, and tabular output columns.

## Non-uniform residue frequencies

By default, CompariMotif uses a uniform residue frequency distribution when it computes
information content. You can override that background model by passing
`residue_frequencies = Dict{Char,Float64}(...)` to `ComparisonOptions`.

```julia
dna_freqs = Dict('A' => 0.3, 'C' => 0.2, 'G' => 0.2, 'T' => 0.3)
weighted = ComparisonOptions(;
    alphabet = DNAAlphabet(),
    residue_frequencies = dna_freqs,
    min_shared_positions = 1,
    normalized_ic_cutoff = 0.0,
)
round(compare("ATG", "[AGT]TG", weighted).match_ic, digits = 3)
# 2.19
```

## Allowed regex symbols and syntax

The supported motif syntax and the corresponding oracle-compatible parser edge cases
are documented on the dedicated `Regex Syntax` page in the docs:

<https://diegozea.github.io/CompariMotif.jl/stable/regex_syntax/>

## Official implementation

The official CompariMotif implementation is distributed as part of SLiMSuite:
<https://github.com/slimsuite/SLiMSuite> (tool path: `tools/comparimotif_V3.py`).

### Scope differences compared to the original CompariMotif

This package implements the paper-defined motif comparison core, but it does not
aim to replicate the full SLiMSuite application surface. In particular:

- no standalone CLI interface or SLiMSuite pipeline integration;
- no raw `.tdt` compatibility/output mode (use `to_column_table` for tabular outputs);
- no `Name*`/`Desc*` metadata fields in API results or fixtures (regex motifs only);
- no XGMML/network export outputs.

#### Implementation notes compared to the paper

This package remains close to the pipeline described in CompariMotif's paper, but a few 
implementation choices differ intentionally:

- No explicit _"enough common amino acids (in any position)"_ prefilter is applied before 
  the sliding-window comparison. Local Julia benchmarking indicated that adding this 
  prefilter increased overall runtime instead of improving performance.
- Exact and exact-subsequence matches do not short-circuit the search. They initialize the 
  current best candidate, but all sliding-window overlaps are still evaluated so that a 
  stronger overlap from another shift or an expanded motif variant can still be selected.
- For partial overlaps between ambiguous residue classes, the wording in the paper can be 
  interpreted as using the lower information content (IC) of the two positions. However, 
  the upstream black-box oracle behaves as if the information content of the union residue 
  class is used instead. This implementation follows the oracle behavior for scoring 
  while retaining the paper’s `Complex` relationship terminology rather than the 
  oracle’s `Ugly` label.
- Ranged repeats and alternations are expanded into concrete motif variants. This 
  expansion is limited by `max_variants` in `ComparisonOptions` (default `10_000`) to 
  prevent pathological combinatorial growth.

### Fixtures and oracle regeneration

Oracle fixtures, i.e. expected results for black-box tests, are committed under 
`data/fixtures/` and tests do not call the CompariMotif code directly. Only normalized 
TSV fixtures are committed rather than the raw `.tdt` output. To regenerate fixtures 
see the `README.md` in `data/fixtures/`.

### Default options parity

Compared against the upstream CompariMotif oracle as a black-box executable
(without reading source code), package defaults match:

- `min_shared_positions = 2` (`minshare=2`)
- `normalized_ic_cutoff = 0.5` (`normcut=0.5`)
- `matchfix = :none` (`matchfix=0`)
- `mismatches = 0`
- `allow_ambiguous_overlap = true` (`overlaps=T`)

### License hygiene

This repository is MIT-licensed. Implementation is derived from the paper and
black-box oracle observations only. GPL CompariMotif source code is not used.

## Citation

If you use this Julia pipeline in scientific work, **please cite CompariMotif's paper**:

- Edwards RJ, Davey NE, Shields DC. *CompariMotif: quick and easy comparisons of
  sequence motifs*. Bioinformatics 24(10):1307-1309 (2008).
  <https://doi.org/10.1093/bioinformatics/btn105>
