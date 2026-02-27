# CompariMotif.jl

Clean-room, unofficial Julia implementation of the motif–motif comparison strategy 
described by Edwards, Davey and Shields (Bioinformatics 24(10):1307–1309, 2008). It 
supports the comparison of protein, DNA and RNA motifs, represented as regular 
expressions.

## API

- `ComparisonOptions(; kwargs...)`
- `compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)::ComparisonResult`
- `compare(motifs::AbstractVector{<:AbstractString}, db::AbstractVector{<:AbstractString}, options::ComparisonOptions)::Matrix{ComparisonResult}`
- `compare(motifs::AbstractVector{<:AbstractString}, options::ComparisonOptions)::Matrix{ComparisonResult}`
- `normalize_motif(motif::AbstractString; alphabet = :protein)::String`
- `to_column_table(result_or_results)::NamedTuple`

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

## Supported motif syntax

- Fixed residues from the selected alphabet:
  - protein (`alphabet=:protein`, default): `ARNDCQEGHILKMFPSTWYV`
  - DNA (`alphabet=:dna`): `ACGT`
  - RNA (`alphabet=:rna`): `ACGU`
- Wildcards: `x`, `X`, and `.`
- Character classes: `[KR]`
- Negated classes: `[^P]` (complement within selected alphabet)
- Anchors: `^` and `$`
- Repeat quantifiers: `{n}`, `{m,n}`, `(n)`, `(m,n)`
- Whitespace is ignored inside motifs.

Unsupported syntax includes general regex alternation/group constructs such as `(A|B)`.

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

### Fixtures and oracle regeneration

Oracle fixtures, i.e. expected results for black-box tests, are committed under 
`data/fixtures/` and tests do not call the CompariMotif code directly. Only normalized 
TSV fixtures are committed rather than the raw `.tdt` output. To regenerate fixtures 
see the `README.md` in `data/fixtures/`.

### License hygiene

This repository is MIT-licensed. Implementation is derived from the paper and
black-box oracle observations only. GPL CompariMotif source code is not used.

## Citation

If you use this Julia pipeline in scientific work, please cite the original
algorithm paper:

- Edwards RJ, Davey NE, Shields DC. *CompariMotif: quick and easy comparisons of
  sequence motifs*. Bioinformatics 24(10):1307-1309 (2008).
  <https://doi.org/10.1093/bioinformatics/btn105>
