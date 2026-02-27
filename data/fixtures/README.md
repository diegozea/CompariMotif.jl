# Fixture Regeneration

This directory stores committed regression fixtures for oracle-based comparisons.

## Files

- `generate_oracle_fixture.jl`: regeneration script (black-box oracle usage only).
- `regression_motifs.tsv`: input motif set used for fixture generation (`#` comments allowed; one motif regex per non-comment line).
- `default_probe_motifs.tsv`: discriminative motif set for default-option parity checks.
- `alternation_probe_motifs.tsv`: motif set with grouping/alternation syntax for oracle parity checks.
- `oracle_regression_normalized.tsv`: normalized, deterministic oracle output used by tests.
- `oracle_default_probe_normalized.tsv`: normalized oracle-default output used by tests.
- `oracle_alternation_probe_normalized.tsv`: normalized oracle output for grouping/alternation tests.

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
```

To regenerate one fixture set only:

```bash
export SLiMSuite_PATH=/path/to/SLiMSuite
julia --project=. data/fixtures/generate_oracle_fixture.jl regression
julia --project=. data/fixtures/generate_oracle_fixture.jl defaults
julia --project=. data/fixtures/generate_oracle_fixture.jl alternation
```

### Oracle Default Probe Fixture

`oracle_default_probe_normalized.tsv` is generated from oracle defaults
(`minshare=2`, `normcut=0.5`, `matchfix=0`, `mismatches=0`, `overlaps=T`)
using `default_probe_motifs.tsv`.

## Test Independence

Package tests do not call the oracle at runtime; they only read the committed normalized fixture.
