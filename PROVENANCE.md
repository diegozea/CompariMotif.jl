# Provenance

## What this repository is

CompariMotif.jl is an independent Julia implementation of the CompariMotif
method described by Edwards, Davey, and Shields (2008). It is an unofficial
package maintained separately from SLiMSuite and should not be described as a
relicensing, continuation, or endorsed port of the original codebase.

## Development approach

The implementation work in this repository was based on two inputs:

- the published method description in Edwards, Davey, and Shields (2008);
- black-box comparison against the behavior of the original `comparimotif_V3.py` script 
  across different motif sets, option settings, and edge cases.

## Sources used

- The CompariMotif paper for the algorithm description, relationship
  terminology, and worked examples.
- The public input/output behavior of the official executable in SLiMSuite,
  treated as an external oracle.
- Oracle-derived fixtures committed under [`data/fixtures/`](data/fixtures/)
  after normalization for deterministic regression testing.

## What was not used

This repository was not produced by copying or translating the original
CompariMotif source code.

Repository code, tests, comments, and documentation were written for this Julia
package rather than imported from the original project. The upstream program may
be used as a black-box executable when regenerating fixtures, but its source is
outside the development inputs for this repository.

## Validation approach

- Unit tests cover parser and scoring primitives.
- Paper examples are reproduced where practical.
- Regression tests compare Julia results with committed fixtures generated from
  black-box runs of the original executable.
- Test runs do not require the upstream tool; the oracle is only needed when
  fixtures are intentionally regenerated.

## Relationship to the original project

The original CompariMotif software distributed with SLiMSuite remains the
official upstream project:
<https://github.com/slimsuite/SLiMSuite>

CompariMotif.jl is a separate Julia implementation of the same published
method. It aims to reproduce the paper-defined comparison behavior within a
Julia library API, not to replace the broader SLiMSuite application surface. No
endorsement by the original authors is implied.

## Motivation for the Julia-native MIT implementation

The motivation for an independent Julia implementation is practical. Much of the Julia 
ecosystem is MIT-licensed, and a native Julia package enables integration with other 
MIT-licensed Julia packages, allowing downstream workflows to call motif comparison in 
process and share Julia data structures directly. Under the original license, routine use 
would likely have had to go through a CLI or subprocess, reducing both 
performance and convenience. The goal is better Julia ecosystem integration, not a 
rebranding or relicensing of the original software.

## Attribution and citation expectations

If you use this package in scientific work, please cite the CompariMotif paper and give 
credit to the original authors:

- Edwards RJ, Davey NE, Shields DC. *CompariMotif: quick and easy comparisons of
  sequence motifs*. Bioinformatics 24(10):1307-1309 (2008).
  <https://doi.org/10.1093/bioinformatics/btn105>