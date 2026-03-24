# Fixture Regeneration

This directory stores committed regression fixtures for oracle-based comparisons.

## Files

- `generate_oracle_fixture.jl`: regeneration script (black-box oracle usage only).
- `regression_motifs.tsv`: input motif set used for fixture generation (`#` comments allowed; one motif regex per non-comment line).
- `default_probe_motifs.tsv`: discriminative motif set for default-option parity checks.
- `boundary_wildcard_probe_motifs.tsv`: motif set targeting clipped boundary-wildcard containments and the corresponding oracle overlap/exact distinctions.
- `matchfix_probe_motifs.tsv`: motif set targeting `matchfix` overlap clipping parity, while retaining oracle-allowed exact-subsequence cases.
- `coreic_probe_motifs.tsv`: motif set targeting `CoreIC` corner cases and anchor/mismatch-sensitive cases; the exact `A..K` / `A..K` wildcard-denominator probe is covered separately by a direct regression test to avoid duplicated motifs that need names for identification.
- `alternation_probe_motifs.tsv`: motif set with grouping/alternation syntax for oracle parity checks.
- `wildcard_alias_divergence_probe_motifs.tsv`: motif set documenting the oracle's top-level grouped wildcard quirk using the stable `.` spelling; package tests separately lock down the intentional `x`/`X` alias equivalence because grouped `x`/`X` probes are unstable under the local Python-3 oracle.
- `exact_prefilter_probe_motifs.tsv`: motif set targeting exact-match, exact-subsequence, and exact-overlap pre-pass parity checks.
- `score_tiebreak_probe_motifs.tsv`: motif set targeting branch/shift competitions that should be resolved by `Score` after `MatchIC` and `MatchPos`.
- `nonuniform_probe_motifs.tsv`: motif set targeting non-uniform background-frequency scoring, including partial ambiguous overlaps that the oracle labels as `Ugly Match`.
- `nonuniform_probe.aafreq.tsv`: strictly positive protein frequency table used for the non-uniform oracle probe.
- `dna_probe_motifs.tsv`: DNA motif set covering the README/docs example plus wildcard, negation, overlap, and anchor corner cases.
- `dna_nonuniform_probe.aafreq.tsv`: DNA frequency table (`A=0.3`, `C=0.2`, `G=0.2`, `T=0.3`) used for the weighted DNA oracle probe.
- `cornercase_probe_motifs.tsv`: motif set targeting whitespace handling and permissive regex corner cases.
- `oracle_regression_normalized.tsv`: normalized, deterministic oracle output used by tests.
- `oracle_default_probe_normalized.tsv`: normalized oracle-default output used by tests.
- `oracle_boundary_wildcard_probe_normalized.tsv`: normalized oracle output for clipped boundary-wildcard containment parity tests.
- `oracle_query_fixed_probe_normalized.tsv`: normalized oracle output for `matchfix=1` (`:query_fixed`) parity tests.
- `oracle_search_fixed_probe_normalized.tsv`: normalized oracle output for `matchfix=2` (`:search_fixed`) parity tests.
- `oracle_both_fixed_probe_normalized.tsv`: normalized oracle output for `matchfix=3` (`:both_fixed`) parity tests.
- `oracle_coreic_probe_normalized.tsv`: normalized oracle output for low-threshold `CoreIC` corner-case parity tests.
- `oracle_coreic_mismatch_probe_normalized.tsv`: normalized oracle output for `mismatches=1` `CoreIC` parity tests.
- `oracle_alternation_probe_normalized.tsv`: normalized oracle output for grouping/alternation tests.
- `oracle_wildcard_alias_divergence_normalized.tsv`: normalized oracle output used to lock down the intentional wildcard-alias divergence documentation.
- `oracle_exact_prefilter_probe_normalized.tsv`: normalized oracle output for exact-match pre-pass parity tests.
- `oracle_score_tiebreak_probe_normalized.tsv`: normalized oracle output for score-based tie-break tests.
- `oracle_nonuniform_probe_normalized.tsv`: normalized oracle output for non-uniform frequency scoring tests.
- `oracle_dna_probe_normalized.tsv`: normalized oracle output for uniform DNA scoring tests.
- `oracle_dna_nonuniform_probe_normalized.tsv`: normalized oracle output for weighted DNA scoring tests.
- `oracle_cornercase_probe_normalized.tsv`: normalized oracle output for whitespace and regex corner-case parity tests.

## Prerequisites

- `python3`
- Local clone of `https://github.com/slimsuite/SLiMSuite`

Configure:

- `SLiMSuite_PATH=/path/to/SLiMSuite`

## Regenerate Fixtures

From repository root:

```bash
export SLiMSuite_PATH=/path/to/SLiMSuite
julia --project=. data/fixtures/generate_oracle_fixture.jl
```

This rewrites all normalized fixture outputs:

```bash
oracle_regression_normalized.tsv
oracle_default_probe_normalized.tsv
oracle_boundary_wildcard_probe_normalized.tsv
oracle_query_fixed_probe_normalized.tsv
oracle_search_fixed_probe_normalized.tsv
oracle_both_fixed_probe_normalized.tsv
oracle_coreic_probe_normalized.tsv
oracle_coreic_mismatch_probe_normalized.tsv
oracle_alternation_probe_normalized.tsv
oracle_wildcard_alias_divergence_normalized.tsv
oracle_exact_prefilter_probe_normalized.tsv
oracle_score_tiebreak_probe_normalized.tsv
oracle_nonuniform_probe_normalized.tsv
oracle_dna_probe_normalized.tsv
oracle_dna_nonuniform_probe_normalized.tsv
oracle_cornercase_probe_normalized.tsv
```

To regenerate one fixture set only:

```bash
export SLiMSuite_PATH=/path/to/SLiMSuite
julia --project=. data/fixtures/generate_oracle_fixture.jl regression
julia --project=. data/fixtures/generate_oracle_fixture.jl defaults
julia --project=. data/fixtures/generate_oracle_fixture.jl boundary_wildcard
julia --project=. data/fixtures/generate_oracle_fixture.jl query_fixed
julia --project=. data/fixtures/generate_oracle_fixture.jl search_fixed
julia --project=. data/fixtures/generate_oracle_fixture.jl both_fixed
julia --project=. data/fixtures/generate_oracle_fixture.jl coreic
julia --project=. data/fixtures/generate_oracle_fixture.jl coreic_mismatch
julia --project=. data/fixtures/generate_oracle_fixture.jl alternation
julia --project=. data/fixtures/generate_oracle_fixture.jl wildcard_alias_divergence
julia --project=. data/fixtures/generate_oracle_fixture.jl exact_prefilter
julia --project=. data/fixtures/generate_oracle_fixture.jl score_tiebreak
julia --project=. data/fixtures/generate_oracle_fixture.jl nonuniform
julia --project=. data/fixtures/generate_oracle_fixture.jl dna
julia --project=. data/fixtures/generate_oracle_fixture.jl dna_nonuniform
julia --project=. data/fixtures/generate_oracle_fixture.jl cornercases
```

### Oracle Default Probe Fixture

`oracle_default_probe_normalized.tsv` is generated from oracle defaults
(`minshare=2`, `normcut=0.5`, `matchfix=0`, `mismatches=0`, `overlaps=T`)
using `default_probe_motifs.tsv`.

### Boundary-Wildcard Probe Fixture

`oracle_boundary_wildcard_probe_normalized.tsv` is generated with
`minshare=1`, `normcut=0`, and `matchfix=0` using
`boundary_wildcard_probe_motifs.tsv`.

### CoreIC Probe Fixtures

`oracle_coreic_probe_normalized.tsv` is generated with `minshare=1`,
`normcut=0`, `mismatches=0`, and `xgmml=F` using
`coreic_probe_motifs.tsv` and keeps the full normalized oracle best-row parity
corpus for that motif set, excluding the exact `A..K` / `A..K` self-match that
the oracle only exposes when duplicate motifs are present in the same input.

`oracle_coreic_mismatch_probe_normalized.tsv` reuses the same motif set with
`mismatches=1` and likewise keeps the full normalized oracle best-row parity
corpus for that motif set, including mismatch-tolerant anchor rows.

The exact `A..K` / `A..K` `CoreIC` values are checked separately in
`test/regression_coreic_oracle.jl` using direct oracle-derived expectations.

## Test Independence

Package tests do not call the oracle at runtime; they only read the committed normalized fixture.
