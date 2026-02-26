# Fixture Regeneration

This directory stores committed regression fixtures for oracle-based comparisons.

## Files

- `regression_motifs.tsv`: input motif set used for fixture generation (`#` comments allowed; one motif regex per non-comment line).
- `oracle_regression_normalized.tsv`: normalized, deterministic oracle output used by tests.
- `generate_oracle_fixture.jl`: regeneration script (black-box oracle usage only).

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

This rewrites `oracle_regression_normalized.tsv`.

## Test Independence

Package tests do not call the oracle at runtime; they only read the committed normalized fixture.
