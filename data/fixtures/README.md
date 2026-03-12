# Fixture Regeneration

This directory stores committed regression fixtures for oracle-based comparisons.

## Files

- `generate_oracle_fixture.jl`: regeneration script (black-box oracle usage only).
- `regression_motifs.tsv`: input motif set used for fixture generation (`#` comments allowed; one motif regex per non-comment line).
- `default_probe_motifs.tsv`: discriminative motif set for default-option parity checks.
- `alternation_probe_motifs.tsv`: motif set with grouping/alternation syntax for oracle parity checks.
- `exact_prefilter_probe_motifs.tsv`: motif set targeting exact-match, exact-subsequence, and exact-overlap pre-pass parity checks.
- `score_tiebreak_probe_motifs.tsv`: motif set targeting branch/shift competitions that should be resolved by `Score` after `MatchIC` and `MatchPos`.
- `nonuniform_probe_motifs.tsv`: motif set targeting non-uniform background-frequency scoring, including partial ambiguous overlaps that the oracle labels as `Ugly Match`.
- `nonuniform_probe.aafreq.tsv`: strictly positive protein frequency table used for the non-uniform oracle probe.
- `oracle_regression_normalized.tsv`: normalized, deterministic oracle output used by tests.
- `oracle_default_probe_normalized.tsv`: normalized oracle-default output used by tests.
- `oracle_alternation_probe_normalized.tsv`: normalized oracle output for grouping/alternation tests.
- `oracle_exact_prefilter_probe_normalized.tsv`: normalized oracle output for exact-match pre-pass parity tests.
- `oracle_score_tiebreak_probe_normalized.tsv`: normalized oracle output for score-based tie-break tests.
- `oracle_nonuniform_probe_normalized.tsv`: normalized oracle output for non-uniform frequency scoring tests.

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
oracle_alternation_probe_normalized.tsv
oracle_exact_prefilter_probe_normalized.tsv
oracle_score_tiebreak_probe_normalized.tsv
oracle_nonuniform_probe_normalized.tsv
```

To regenerate one fixture set only:

```bash
export SLiMSuite_PATH=/path/to/SLiMSuite
julia --project=. data/fixtures/generate_oracle_fixture.jl regression
julia --project=. data/fixtures/generate_oracle_fixture.jl defaults
julia --project=. data/fixtures/generate_oracle_fixture.jl alternation
julia --project=. data/fixtures/generate_oracle_fixture.jl exact_prefilter
julia --project=. data/fixtures/generate_oracle_fixture.jl score_tiebreak
julia --project=. data/fixtures/generate_oracle_fixture.jl nonuniform
```

### Oracle Default Probe Fixture

`oracle_default_probe_normalized.tsv` is generated from oracle defaults
(`minshare=2`, `normcut=0.5`, `matchfix=0`, `mismatches=0`, `overlaps=T`)
using `default_probe_motifs.tsv`.

## Test Independence

Package tests do not call the oracle at runtime; they only read the committed normalized fixture.
