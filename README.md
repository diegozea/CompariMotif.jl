# CompariMotif.jl

Clean-room, unofficial Julia implementation of the motif–motif comparison strategy 
described by Edwards, Davey and Shields (Bioinformatics 24(10):1307–1309, 2008). It 
supports the comparison of protein and DNA motifs, represented as regular 
expressions.

## API

- `ComparisonOptions(; kwargs...)`
- `compare(a::AbstractString, b::AbstractString, options::ComparisonOptions)::ComparisonResult`
- `compare(motifs::AbstractVector{<:AbstractString}, db::AbstractVector{<:AbstractString}, options::ComparisonOptions)::Matrix{ComparisonResult}`
- `compare(motifs::AbstractVector{<:AbstractString}, options::ComparisonOptions)::Matrix{ComparisonResult}`
- `normalize_motif(motif::AbstractString; alphabet = :protein)::String`
- `write_results_tsv(path, motifs, db, results)`

## Minimal example

```julia
using CompariMotif

motifs = ["RKLI", "R[KR]L[IV]", "[KR]xLx[FYLIMVP]", "RxLE"]
options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0)
results = compare(motifs, options)

results[3, 4]  # single pair summary
write_results_tsv("comparimotif_results.tsv", motifs, motifs, results)  # save full matrix
```

## Official implementation

The official CompariMotif implementation is distributed as part of SLiMSuite:
<https://github.com/slimsuite/SLiMSuite> (tool path: `tools/comparimotif_V3.py`).

### Scope differences compared to the original CompariMotif

This package implements the paper-defined motif comparison core, but it does not
aim to replicate the full SLiMSuite application surface. In particular:

- no standalone CLI interface or SLiMSuite pipeline integration;
- no raw `.tdt` compatibility/output mode (deterministic TSV writer is provided);
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
