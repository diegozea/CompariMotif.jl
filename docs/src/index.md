# CompariMotif

CompariMotif.jl is a clean-room Julia implementation of the motif comparison workflow
described in Edwards et al. (2008).

## Quick Start

```jldoctest
julia> using CompariMotif

julia> motifs = ["RKLI", "R[KR]L[IV]", "[KR]xLx[FYLIMVP]", "RxLE"];

julia> results = compare(motifs; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> (results[3, 4].matched, results[3, 4].query_relationship, results[3, 4].matched_positions)
(true, "Degenerate Parent", 2)

julia> outfile = joinpath(mktempdir(), "comparimotif_results.tsv");

julia> write_results_tsv(outfile, motifs, motifs, results);

julia> isfile(outfile)
true
```

## Public API

```@autodocs
Modules = [CompariMotif]
Private = false
```
