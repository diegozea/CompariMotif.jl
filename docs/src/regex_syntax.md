```@meta
CurrentModule = Main.CompariMotif
```

# Regex Syntax

Motif parsing supports a controlled regex-like subset.

- Fixed residues from the selected alphabet:
  - protein (`alphabet=ProteinAlphabet()`, default): `ARNDCQEGHILKMFPSTWYV`
  - DNA (`alphabet=DNAAlphabet()`): `ACGT`
  - RNA (`alphabet=RNAAlphabet()`): `ACGU`
- Wildcards:
  - `x`, `X`, and `.` are equivalent and mean "any residue in the selected alphabet".
- Character classes:
  - `[KR]` includes listed residues.
  - `[^P]` is negation within the selected alphabet only.
- Anchors:
  - `^` and `$` indicate N- and C-terminus for protein motifs, or 5' and 3' ends for nucleic acid motifs.
- Repeat quantifiers:
  - `{n}`, `{m,n}` on residues, classes, anchors, and grouped alternatives.
- Grouping and alternation:
  - `(...)` for grouping and `|` for alternatives, for example `A(K|Q)LI`.
- Whitespace:
  - leading and trailing whitespace is trimmed;
  - the first internal whitespace ends motif parsing, matching the upstream oracle.

## Edge Cases and Oracle Behavior

For a handful of parser edge cases, this package intentionally matches the observed
black-box behavior of the upstream CompariMotif oracle rather than enforcing a
stricter regex grammar.

- Internal whitespace truncates parsing at the first space, so `A C` behaves as `A`.
- Some permissive exact quantifiers are accepted and normalized the same way as the
  oracle, for example `A{0}` and `A{-1}` both behave as `A`.
- Grouped exact quantifiers are supported to mirror the oracle split behavior, for
  example `(A|C){2}` expands to `AA` and `CC`.
- Some malformed constructs are treated as non-retained motifs rather than hard parse
  failures, such as `A{2,1}` producing zero retained variants.
- Truly malformed syntax is still rejected when the oracle also rejects it, for
  example an unclosed character class like `A[Q`.
